#!/usr/bin/env python3
"""run_protein.py - one command, one protein.

    python run_protein.py STRUCTURE

Give it a protein structure (a PDB/mmCIF file, or a 4-character PDB ID) and it
produces, in the output folder:

    NAME.graphml        the residue network
    NAME.qmod           the Classiq quantum-walk circuit
    NAME_hump.png       the graph: signal transport vs noise (look for the hump)
    NAME_hump.csv       the numbers behind that graph
    NAME_execute.txt    how to run the .qmod on Classiq

The hump graph is the point: it sweeps the dephasing rate and shows how much
signal travels from a source residue to a distant one. A peak at intermediate
noise is the environment-assisted effect. No labels or databases needed.

Options you may want: --source to set the starting residue (default: the
most-connected one), --site-energy / --scale for the energy model, and the usual
network options (--chains, --method, --cutoff). --no-qmod skips the circuit if
you only want the graph.
"""
import argparse, os, sys
import numpy as np, networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rin_builder as rb
import walk_core as wc

EXEC_SCAFFOLD = '''\
# Run this .qmod on Classiq after classiq.authenticate().
# Noise-free (quantum end):
#   qprog = synthesize(create_model(main)); execute(qprog).result()
# Environment-assisted (the hump): run on SIMULATOR_DENSITY_MATRIX with a
# dephasing channel whose rate is your gamma dial. Confirm the exact
# ClassiqSimulatorNoiseSpecification fields on the live "custom noise models"
# docs page, then sweep the rate and check it matches the hump graph here.
'''


def structure_to_graph(inp, outdir, model, chains, method, cutoff, min_seq_sep, weighted):
    residues = rb.collect_residues(rb.load_structure(inp, tmpdir=outdir), model, chains)
    if not residues:
        raise ValueError("no amino-acid residues found.")
    parts = (rb.build_by_heavy_atoms(residues, cutoff, min_seq_sep) if method == "heavy"
             else rb.build_by_representative(residues, method, cutoff, min_seq_sep))
    return rb.assemble_graph(*parts, weighted)


def write_qmod(terms, n_qubits, source_index, evolution_time, order, repetitions, out_path):
    """Emit the Qmod (offline, no login). Source is amplitude-encoded with qubit
    0 as the least significant bit, matching Classiq's convention."""
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
    ap = argparse.ArgumentParser(description="One protein: structure -> Qmod + hump graph.")
    ap.add_argument("input", help="PDB/mmCIF file, or a 4-character PDB ID.")
    ap.add_argument("--source", default=None, help="Start residue id, e.g. A:151 (default: most connected).")
    ap.add_argument("--target", default=None, help="Transport target (default: farthest residue).")
    # network options
    ap.add_argument("--chains", default=None, help="Chains to include (default all).")
    ap.add_argument("--model", type=int, default=0, help="Model index for NMR files.")
    ap.add_argument("--method", choices=["cb", "ca", "heavy"], default="cb")
    ap.add_argument("--cutoff", type=float, default=None)
    ap.add_argument("--min-seq-sep", type=int, default=0)
    ap.add_argument("--weighted", action="store_true")
    # physics options
    ap.add_argument("--site-energy", choices=["none", "hydropathy", "degree", "random"], default="hydropathy")
    ap.add_argument("--scale", type=float, default=3.0, help="Site-energy disorder strength.")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--gammas", default="0,0.05,0.1,0.2,0.4,0.8,1.5,3,6,12,25,50")
    ap.add_argument("--tmax", type=float, default=30.0)
    ap.add_argument("--ntime", type=int, default=200)
    # qmod options
    ap.add_argument("--time", type=float, default=1.0, help="Evolution time in the circuit.")
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--repetitions", type=int, default=4)
    ap.add_argument("--max-residues", type=int, default=200, help="Skip the Qmod above this size unless --force.")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-qmod", action="store_true", help="Only make the hump graph.")
    ap.add_argument("--outdir", default="output")
    ap.add_argument("--prefix", default=None)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    prefix = a.prefix or (os.path.splitext(os.path.basename(a.input))[0] or "protein")
    cutoff = a.cutoff if a.cutoff is not None else (5.0 if a.method == "heavy" else 8.0)
    chains = {c.strip() for c in a.chains.split(",")} if a.chains else None

    # 1) structure -> network
    print(f"[1/3] structure -> network: {a.input}")
    try:
        G = structure_to_graph(a.input, a.outdir, a.model, chains, a.method, cutoff, a.min_seq_sep, a.weighted)
    except Exception as e:
        sys.exit(f"\n[error] could not build the network from '{a.input}': {e}\n"
                 f"If you gave a PDB ID, check the ID and your internet connection.")
    graphml = os.path.join(a.outdir, prefix + ".graphml")
    nx.write_graphml(G, graphml)
    _, nodes, A, resnames = wc.load_network(graphml)
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    n_qubits = max(1, int(np.ceil(np.log2(N))))
    print(f"      residues={N}  chains={sorted({n.split(':')[0] for n in nodes})}  -> qubits={n_qubits}")

    # pick source and target
    if a.source and a.source not in idx:
        sys.exit(f"[error] source '{a.source}' not in network, e.g. {nodes[:5]}")
    source_id = a.source or max(dict(G.degree()), key=dict(G.degree()).get)
    if a.target and a.target in idx:
        target_id = a.target
    else:
        target_id = max(nx.single_source_shortest_path_length(G, source_id).items(),
                        key=lambda kv: kv[1])[0]
    print(f"      source {source_id}"
          + ("" if a.source else " (auto: most connected)")
          + f"  -> target {target_id}"
          + ("" if a.target else " (auto: farthest)"))

    # 2) hump graph (classical, exact, any size)
    print(f"[2/3] transport vs noise ({a.site_energy} site energies, scale {a.scale})")
    H = wc.build_hamiltonian(A, resnames, a.site_energy, a.scale, a.seed)
    gammas = [float(x) for x in a.gammas.split(",")]
    tlist = np.linspace(0, a.tmax, a.ntime)
    transport = [wc.visiting_scores(wc.run_walk(H, idx[source_id], g, tlist), tlist)[idx[target_id]]
                 for g in gammas]
    peak = int(np.argmax(transport))

    with open(os.path.join(a.outdir, prefix + "_hump.csv"), "w") as f:
        f.write("gamma,transport\n")
        f.writelines(f"{g},{v:.6f}\n" for g, v in zip(gammas, transport))
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 5))
        plt.plot(gammas, transport, "o-", lw=2)
        plt.axvline(gammas[peak], color="crimson", ls="--", alpha=0.7, label=f"peak at gamma={gammas[peak]:g}")
        plt.xlabel("dephasing rate gamma  (0 = fully quantum, large = classical)")
        plt.ylabel(f"signal reaching {target_id}")
        plt.title(f"{prefix}: transport vs noise\nsource {source_id} -> target {target_id}")
        plt.legend(); plt.tight_layout()
        plt.savefig(os.path.join(a.outdir, prefix + "_hump.png"), dpi=140); plt.close()
    except Exception as e:
        print(f"[warn] plot skipped: {e}")
    verdict = ("HUMP: transport peaks at intermediate noise (environment-assisted)."
               if 0 < peak < len(gammas) - 1 else
               "No hump: best fully quantum (gamma=0)." if peak == 0 else
               "No hump: best at the classical end.")
    print(f"      {verdict}")

    # 3) qmod
    if a.no_qmod:
        print("[3/3] Qmod skipped (--no-qmod)")
    elif N > a.max_residues and not a.force:
        print(f"[3/3] Qmod skipped: {N} residues > --max-residues={a.max_residues} "
              f"(the circuit would be impractically large). Use --chains, a single "
              f"domain, or --force. The hump graph above is unaffected.")
    else:
        print(f"[3/3] Hamiltonian -> Pauli -> Qmod")
        Hp, nq = wc.pad_to_power_of_two(H)
        terms, n = wc.pauli_decompose(Hp)
        err = float(np.max(np.abs(wc.pauli_reconstruct(terms, n) - Hp)))
        out = os.path.join(a.outdir, prefix + ".qmod")
        qmod_file, n_terms = write_qmod(terms, nq, idx[source_id], a.time, a.order, a.repetitions, out)
        open(os.path.join(a.outdir, prefix + "_execute.txt"), "w").write(EXEC_SCAFFOLD)
        print(f"      {n_terms} Pauli terms (decomp error {err:.1e})  ->  {qmod_file} "
              f"({os.path.getsize(qmod_file)/1024:.0f} KB)")

    print(f"\nDone. Outputs in {a.outdir}/")
    print(f"  hump graph : {os.path.join(a.outdir, prefix + '_hump.png')}")
    print(f"  network    : {graphml}")
    if not a.no_qmod and (N <= a.max_residues or a.force):
        print(f"  qmod       : {os.path.join(a.outdir, prefix + '.qmod')}")


if __name__ == "__main__":
    main()
