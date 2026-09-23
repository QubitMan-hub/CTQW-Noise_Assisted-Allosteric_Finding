#!/usr/bin/env python3
"""One command: any protein structure -> a Classiq Qmod walk file.
Chains rin_builder (structure -> network) and walk_core + qmod_walk (network ->
Qmod), with PDB-ID or file input, chain selection, a default source, and a size
guard against uselessly large Qmods.

    python protein_to_qmod.py 1A8O.pdb --source A:151
    python protein_to_qmod.py 6LU7 --chains A
"""
import argparse, os, sys
import numpy as np, networkx as nx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import rin_builder as rb, walk_core as wc, qmod_walk as qw
except ImportError as e:
    sys.exit(f"[error] keep rin_builder.py, walk_core.py, qmod_walk.py beside this file.\n  {e}")


def structure_to_graph(inp, outdir, model, chains, method, cutoff, min_seq_sep, weighted):
    residues = rb.collect_residues(rb.load_structure(inp, tmpdir=outdir), model, chains)
    if not residues:
        raise ValueError("no amino-acid residues found.")
    parts = (rb.build_by_heavy_atoms(residues, cutoff, min_seq_sep) if method == "heavy"
             else rb.build_by_representative(residues, method, cutoff, min_seq_sep))
    return rb.assemble_graph(*parts, weighted)


def main():
    ap = argparse.ArgumentParser(description="Any protein structure -> Qmod walk.")
    ap.add_argument("input", help="PDB/mmCIF file or 4-character PDB ID.")
    ap.add_argument("--source", default=None, help="Source residue id (default: highest-degree).")
    ap.add_argument("--chains", default=None)
    ap.add_argument("--model", type=int, default=0)
    ap.add_argument("--method", choices=["cb", "ca", "heavy"], default="cb")
    ap.add_argument("--cutoff", type=float, default=None)
    ap.add_argument("--min-seq-sep", type=int, default=0)
    ap.add_argument("--weighted", action="store_true")
    ap.add_argument("--site-energy", choices=["none", "hydropathy", "degree", "random"], default="hydropathy")
    ap.add_argument("--scale", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--time", type=float, default=1.0)
    ap.add_argument("--order", type=int, default=2)
    ap.add_argument("--repetitions", type=int, default=4)
    ap.add_argument("--max-residues", type=int, default=200)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--outdir", default="pipeline_output")
    ap.add_argument("--prefix", default=None)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    prefix = a.prefix or (os.path.splitext(os.path.basename(a.input))[0] or "protein")
    cutoff = a.cutoff if a.cutoff is not None else (5.0 if a.method == "heavy" else 8.0)
    chains = {c.strip() for c in a.chains.split(",")} if a.chains else None

    print(f"[1/4] structure -> network: {a.input}")
    try:
        G = structure_to_graph(a.input, a.outdir, a.model, chains, a.method, cutoff, a.min_seq_sep, a.weighted)
    except Exception as e:
        sys.exit(f"\n[error] could not build the network from '{a.input}': {e}\n"
                 f"If you gave a PDB ID, check the ID and your internet connection.")

    graphml = os.path.join(a.outdir, prefix + ".graphml")
    nx.write_graphml(G, graphml)
    _, nodes, A, resnames = wc.load_network(graphml)
    N = len(nodes)
    n_qubits = max(1, int(np.ceil(np.log2(N))))
    print(f"      residues={N}  chains={sorted({n.split(':')[0] for n in nodes})}  -> qubits={n_qubits}")

    if N > a.max_residues and not a.force:
        sys.exit(f"\n[stop] {N} residues exceeds --max-residues={a.max_residues}; the Qmod would be "
                 f"impractically large.\nUse --chains to pick one chain, a single-domain file, --force "
                 f"to override, or the classical route (score_allostery.py / run_benchmark.py, no size limit).")

    idx = {n: i for i, n in enumerate(nodes)}
    if a.source and a.source not in idx:
        sys.exit(f"[error] source '{a.source}' not in network, e.g. {nodes[:5]}")
    source_id = a.source or max(dict(G.degree()), key=dict(G.degree()).get)
    print(f"[2/4] source residue: {source_id}" + ("" if a.source else "  (auto highest-degree; pass --source)"))

    print(f"[3/4] Hamiltonian ({a.site_energy}) -> Pauli decomposition")
    H = wc.build_hamiltonian(A, resnames, a.site_energy, a.scale, a.seed)
    Hp, nq = wc.pad_to_power_of_two(H)
    terms, n = wc.pauli_decompose(Hp)
    err = float(np.max(np.abs(wc.pauli_reconstruct(terms, n) - Hp)))
    print(f"      pauli_terms={len(terms)}  decomposition_error={err:.1e} ({'OK' if err < 1e-6 else 'WARNING'})")

    print("[4/4] writing Qmod")
    out = os.path.join(a.outdir, prefix + ".qmod")
    qmod_file, n_terms = qw.write_qmod(terms, nq, idx[source_id], a.time, a.order, a.repetitions, out)
    open(os.path.join(a.outdir, prefix + "_execute.py.txt"), "w").write(qw.EXEC_SCAFFOLD)
    size_kb = os.path.getsize(qmod_file) / 1024
    print(f"\nDone. {N} residues -> {nq} qubits -> {n_terms} Pauli terms.")
    print(f"Qmod: {qmod_file} ({size_kb:.0f} KB)  |  Network: {graphml}")
    if size_kb > 5000:
        print("[note] large Qmod; synthesis may be slow. Consider one chain/domain.")


if __name__ == "__main__":
    main()
