#!/usr/bin/env python3
"""Network -> Classiq Qmod quantum walk file (.qmod written offline, no login).

Amplitude-encodes residues into ceil(log2 N) qubits, decomposes H into Pauli
terms, and emits a suzuki_trotter circuit that prepares the source and evolves
it (exp(-iHt)). Noise-free; add dephasing at execution time on the density-matrix
simulator (see the printed scaffold). Small proteins only: an arbitrary H gives
up to N^2 Pauli terms.

    python qmod_walk.py NET.graphml --source A:151
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import walk_core as wc

EXEC_SCAFFOLD = '''\
# Run on Classiq after classiq.authenticate().
# Noise-free (quantum end): qprog = synthesize(create_model(main)); execute(qprog).result()
# Environment-assisted (the paper): run on SIMULATOR_DENSITY_MATRIX with a
# dephasing channel whose rate is your gamma dial. Confirm the exact
# ClassiqSimulatorNoiseSpecification fields on the live "custom noise models"
# docs page, then sweep the rate and check it reproduces classical_qsw.py.
'''


def write_qmod(terms, n_qubits, source_index, evolution_time, order, repetitions, out_path):
    from classiq import (Output, QArray, QBit, qfunc, allocate, X,
                         suzuki_trotter, create_model, Pauli, PauliTerm)
    pmap = {"I": Pauli.I, "X": Pauli.X, "Y": Pauli.Y, "Z": Pauli.Z}
    ham = [PauliTerm(pauli=[pmap[c] for c in s], coefficient=v) for s, v in terms]
    set_bits = [q for q in range(n_qubits) if (source_index >> q) & 1]

    @qfunc
    def main(q: Output[QArray[QBit]]):
        allocate(n_qubits, q)
        for b in set_bits:
            X(q[b])
        suzuki_trotter(ham, evolution_coefficient=evolution_time,
                       order=order, repetitions=repetitions, qbv=q)

    base = os.path.splitext(out_path)[0]
    create_model(main, out_file=base)
    return base + ".qmod", len(ham)


def main():
    ap = argparse.ArgumentParser(description="Residue network -> Classiq Qmod walk.")
    ap.add_argument("network")
    ap.add_argument("--source", required=True, help="Source residue id, e.g. A:151")
    ap.add_argument("--site-energy", choices=["none", "hydropathy", "degree", "random"], default="hydropathy")
    ap.add_argument("--scale", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--time", type=float, default=1.0)
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--repetitions", type=int, default=4)
    ap.add_argument("--outdir", default="qmod_output")
    ap.add_argument("--prefix", default=None)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    prefix = a.prefix or os.path.splitext(os.path.basename(a.network))[0] + "_walk"
    G, nodes, A, resnames = wc.load_network(a.network)
    idx = {n: i for i, n in enumerate(nodes)}
    if a.source not in idx:
        sys.exit(f"[error] source '{a.source}' not in network, e.g. {nodes[:5]}")

    H = wc.build_hamiltonian(A, resnames, a.site_energy, a.scale, a.seed)
    Hp, n_qubits = wc.pad_to_power_of_two(H)
    terms, n = wc.pauli_decompose(Hp)
    err = float(np.max(np.abs(wc.pauli_reconstruct(terms, n) - Hp)))
    print(f"residues={len(nodes)}  qubits={n_qubits}  pauli_terms={len(terms)}  "
          f"decomp_error={err:.1e} ({'OK' if err < 1e-6 else 'WARNING'})")

    out = os.path.join(a.outdir, prefix + ".qmod")
    qmod_file, n_terms = write_qmod(terms, n_qubits, idx[a.source], a.time, a.order, a.repetitions, out)
    open(os.path.join(a.outdir, prefix + "_execute.py.txt"), "w").write(EXEC_SCAFFOLD)
    print(f"wrote {qmod_file} ({n_terms} terms; source {a.source} at basis index {idx[a.source]})")


if __name__ == "__main__":
    main()
