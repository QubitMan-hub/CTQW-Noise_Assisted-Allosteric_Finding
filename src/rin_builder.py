#!/usr/bin/env python3
"""
rin_builder.py
Convert a protein structure into a residue interaction network (a graph).

Input : a PDB/mmCIF file, OR a 4-character PDB ID (fetched from RCSB).
Output: the graph saved as GraphML + a sparse adjacency matrix + node/edge
        tables + a text summary, and (optionally) a quick picture.

Nodes are amino-acid residues. An edge joins two residues that are in contact,
where "contact" is defined one of three ways (see --method).

Design goal: run on real-world structures without hand-holding. It handles
multi-model files (NMR), multiple chains, waters and other hetero groups,
alternate conformations, missing atoms, insertion codes, and modified
residues such as selenomethionine (MSE -> MET).

Author-agnostic, pure Python, works the same on Windows/macOS/Linux.
"""

import argparse
import os
import sys
import warnings

import numpy as np
import networkx as nx
from scipy import sparse
from scipy.spatial import cKDTree

from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter("ignore", PDBConstructionWarning)

# Modified residues we treat as their standard parent so they aren't dropped.
MODIFIED_TO_STANDARD = {
    "MSE": "MET",  # selenomethionine
    "SEP": "SER", "TPO": "THR", "PTR": "TYR",  # phosphorylated
    "HYP": "PRO", "PCA": "GLN", "CSO": "CYS",
    "MLY": "LYS", "M3L": "LYS", "KCX": "LYS",
    "CME": "CYS", "CSD": "CYS", "OCS": "CYS", "FME": "MET",
}


def load_structure(inp, tmpdir="."):
    """Return a Bio.PDB Structure from a file path or a 4-char PDB ID."""
    # If it's an existing file, parse by extension.
    if os.path.isfile(inp):
        ext = os.path.splitext(inp)[1].lower()
        if ext in (".cif", ".mmcif"):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        return parser.get_structure("structure", inp)

    # Otherwise treat it as a PDB ID and fetch it.
    token = inp.strip()
    if len(token) == 4 and token.isalnum():
        from Bio.PDB import PDBList
        pdbl = PDBList()
        path = pdbl.retrieve_pdb_file(token, pdir=tmpdir, file_format="pdb")
        if not os.path.isfile(path):  # some IDs only exist as mmCIF
            path = pdbl.retrieve_pdb_file(token, pdir=tmpdir, file_format="mmCif")
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        return parser.get_structure(token, path)

    raise FileNotFoundError(
        f"'{inp}' is neither a readable file nor a 4-character PDB ID."
    )


def standard_resname(residue):
    """Return the standard 3-letter name if this residue is (or maps to) an
    amino acid, else None."""
    name = residue.get_resname().strip().upper()
    if name in MODIFIED_TO_STANDARD:
        return MODIFIED_TO_STANDARD[name]
    # is_aa with standard=True catches the 20 canonical residues.
    if is_aa(name, standard=True):
        return name
    return None


def pick_representative_atom(residue, method):
    """Return (name, coord) for the atom that represents this residue, or None.
    For 'cb' we use CB, falling back to CA (glycine has no CB). For 'ca' we use
    CA. Alternate conformations: Biopython gives the first/selected altloc."""
    def coord_of(atom_name):
        if atom_name in residue:
            return residue[atom_name].get_coord()
        return None

    if method == "ca":
        c = coord_of("CA")
        return ("CA", c) if c is not None else None
    # method == 'cb'
    c = coord_of("CB")
    if c is not None:
        return ("CB", c)
    c = coord_of("CA")  # glycine and any residue missing CB
    return ("CA", c) if c is not None else None


def collect_residues(structure, model_index, chains_wanted):
    """Walk the chosen model and return an ordered list of accepted residues,
    each as a dict of identity + the Residue object."""
    models = list(structure.get_models())
    if not models:
        raise ValueError("Structure contains no models.")
    if model_index >= len(models):
        raise ValueError(
            f"Requested model {model_index} but only {len(models)} present."
        )
    model = models[model_index]

    residues = []
    for chain in model:
        cid = chain.get_id().strip() or "_"
        if chains_wanted and cid not in chains_wanted:
            continue
        for res in chain:
            hetflag, resseq, icode = res.get_id()
            std = standard_resname(res)
            if std is None:
                continue  # water, ligand, ion, nucleic acid, etc.
            residues.append({
                "chain": cid,
                "resseq": resseq,
                "icode": (icode or "").strip(),
                "resname": std,
                "orig_resname": res.get_resname().strip(),
                "residue": res,
            })
    return residues


def node_label(info):
    """Stable, human-readable node id like 'A:57' or 'A:57B' (insertion code)."""
    return f"{info['chain']}:{info['resseq']}{info['icode']}"


def build_by_representative(residues, method, cutoff, min_seq_sep):
    """Contacts via a single representative atom per residue (CB/CA)."""
    labels, coords, kept = [], [], []
    for info in residues:
        rep = pick_representative_atom(info["residue"], method)
        if rep is None or rep[1] is None:
            continue  # residue has neither CB nor CA; skip it
        labels.append(node_label(info))
        coords.append(rep[1])
        kept.append(info)

    coords = np.asarray(coords, dtype=float)
    edges = []
    if len(coords) >= 2:
        tree = cKDTree(coords)
        pairs = tree.query_pairs(r=cutoff)  # set of (i, j), i < j
        for i, j in pairs:
            if _seq_far_enough(kept[i], kept[j], min_seq_sep):
                d = float(np.linalg.norm(coords[i] - coords[j]))
                edges.append((labels[i], labels[j], d))
    return labels, kept, coords, edges


def build_by_heavy_atoms(residues, cutoff, min_seq_sep):
    """Contacts via the minimum distance between any pair of heavy (non-H)
    atoms of two residues. More physical, a little slower."""
    atom_list, owner = [], []  # owner[k] = index into `residues` for atom k
    labels = [node_label(info) for info in residues]
    for idx, info in enumerate(residues):
        for atom in info["residue"]:
            el = (atom.element or atom.get_name()[0]).strip().upper()
            if el == "H" or atom.get_name().startswith("H"):
                continue  # skip hydrogens
            atom_list.append(atom)
            owner.append(idx)

    edges = {}
    if len(atom_list) >= 2:
        ns = NeighborSearch(atom_list)
        close = ns.search_all(cutoff, level="A")  # all atom pairs within cutoff
        # Map each atom object to the residue it belongs to (fast lookup).
        atom_to_owner = {id(a): owner[k] for k, a in enumerate(atom_list)}
        for a, b in close:
            ia, ib = atom_to_owner[id(a)], atom_to_owner[id(b)]
            if ia == ib:
                continue
            if not _seq_far_enough(residues[ia], residues[ib], min_seq_sep):
                continue
            key = (min(ia, ib), max(ia, ib))
            d = float(np.linalg.norm(a.get_coord() - b.get_coord()))
            if key not in edges or d < edges[key]:
                edges[key] = d

    coords = []
    for info in residues:
        rep = pick_representative_atom(info["residue"], "ca")
        coords.append(rep[1] if rep and rep[1] is not None else (np.nan,)*3)
    coords = np.asarray(coords, dtype=float)

    edge_list = [(labels[i], labels[j], d) for (i, j), d in edges.items()]
    return labels, residues, coords, edge_list


def _seq_far_enough(info_a, info_b, min_seq_sep):
    """True unless the two residues are on the same chain and closer than
    min_seq_sep apart in sequence. Different chains always pass."""
    if min_seq_sep <= 0:
        return True
    if info_a["chain"] != info_b["chain"]:
        return True
    return abs(info_a["resseq"] - info_b["resseq"]) >= min_seq_sep


def assemble_graph(labels, kept, coords, edges, weighted):
    """Build a NetworkX graph with node attributes and edge distances."""
    G = nx.Graph()
    for k, info in enumerate(kept):
        x, y, z = (float(coords[k][0]), float(coords[k][1]), float(coords[k][2])) \
            if coords is not None and not np.any(np.isnan(coords[k])) else (0.0, 0.0, 0.0)
        G.add_node(
            labels[k],
            chain=info["chain"],
            resseq=int(info["resseq"]),
            icode=info["icode"],
            resname=info["resname"],
            orig_resname=info["orig_resname"],
            x=x, y=y, z=z,
        )
    for u, v, d in edges:
        w = (1.0 / d) if (weighted and d > 0) else 1.0
        G.add_edge(u, v, distance=round(d, 3), weight=round(w, 6))
    return G


def write_outputs(G, outdir, prefix, method, cutoff, weighted, want_plot):
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, prefix)

    # 1) GraphML: the full graph with all attributes, reloadable anywhere.
    nx.write_graphml(G, base + ".graphml")

    # 2) Ordered node table (index -> residue identity) so the adjacency
    #    matrix rows/columns are unambiguous.
    nodes = list(G.nodes(data=True))
    with open(base + "_nodes.csv", "w") as fh:
        fh.write("index,node,chain,resseq,icode,resname,x,y,z\n")
        for i, (n, a) in enumerate(nodes):
            fh.write(f"{i},{n},{a['chain']},{a['resseq']},{a['icode']},"
                     f"{a['resname']},{a['x']:.3f},{a['y']:.3f},{a['z']:.3f}\n")

    # 3) Edge table.
    with open(base + "_edges.csv", "w") as fh:
        fh.write("source,target,distance,weight\n")
        for u, v, a in G.edges(data=True):
            fh.write(f"{u},{v},{a['distance']},{a['weight']}\n")

    # 4) Sparse adjacency matrix (weighted if requested, else binary).
    order = [n for n, _ in nodes]
    A = nx.to_scipy_sparse_array(
        G, nodelist=order, weight=("weight" if weighted else None), format="csr"
    )
    sparse.save_npz(base + "_adjacency.npz", A)

    # 5) Human-readable summary.
    n, m = G.number_of_nodes(), G.number_of_edges()
    dens = nx.density(G) if n > 1 else 0.0
    ncomp = nx.number_connected_components(G) if n else 0
    degs = [d for _, d in G.degree()]
    with open(base + "_summary.txt", "w") as fh:
        fh.write("Residue interaction network summary\n")
        fh.write("-----------------------------------\n")
        fh.write(f"contact method     : {method}\n")
        fh.write(f"distance cutoff (A) : {cutoff}\n")
        fh.write(f"edge weighting      : {'1/distance' if weighted else 'binary'}\n")
        fh.write(f"residues (nodes)    : {n}\n")
        fh.write(f"contacts (edges)    : {m}\n")
        fh.write(f"graph density       : {dens:.4f}\n")
        fh.write(f"connected comps     : {ncomp}\n")
        if degs:
            fh.write(f"mean/min/max degree : "
                     f"{np.mean(degs):.2f} / {min(degs)} / {max(degs)}\n")
        chains = sorted({a['chain'] for _, a in nodes})
        fh.write(f"chains included     : {', '.join(chains)}\n")

    # 6) Optional picture.
    if want_plot and n:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            pos = nx.spring_layout(G, seed=1, k=None)
            deg = np.array(degs, dtype=float)
            plt.figure(figsize=(9, 9))
            nx.draw_networkx_edges(G, pos, alpha=0.2, width=0.5)
            nx.draw_networkx_nodes(G, pos, node_size=25 + 4 * deg,
                                   node_color=deg, cmap="viridis")
            plt.axis("off")
            plt.title(f"{prefix}: {n} residues, {m} contacts")
            plt.tight_layout()
            plt.savefig(base + "_network.png", dpi=140)
            plt.close()
        except Exception as e:
            print(f"[warn] plot skipped: {e}", file=sys.stderr)

    return base


def main():
    ap = argparse.ArgumentParser(
        description="Convert a protein structure into a residue interaction "
                    "network (graph).")
    ap.add_argument("input", help="PDB/mmCIF file path, or a 4-character PDB ID.")
    ap.add_argument("--method", choices=["cb", "ca", "heavy"], default="cb",
                    help="Contact definition. cb: CB-CB (CA for Gly). "
                         "ca: CA-CA. heavy: min distance between heavy atoms. "
                         "Default: cb.")
    ap.add_argument("--cutoff", type=float, default=None,
                    help="Contact distance cutoff in angstrom. "
                         "Default 8.0 for cb/ca, 5.0 for heavy.")
    ap.add_argument("--chains", default=None,
                    help="Comma-separated chain IDs to include (default: all).")
    ap.add_argument("--model", type=int, default=0,
                    help="Model index for multi-model (NMR) files. Default: 0.")
    ap.add_argument("--min-seq-sep", type=int, default=0,
                    help="Drop same-chain contacts closer than this in sequence "
                         "(e.g. 2 removes trivial backbone neighbours). Default 0.")
    ap.add_argument("--weighted", action="store_true",
                    help="Write a 1/distance weighted adjacency instead of binary.")
    ap.add_argument("--outdir", default="rin_output", help="Output directory.")
    ap.add_argument("--prefix", default=None,
                    help="Output filename prefix. Default: derived from input.")
    ap.add_argument("--no-plot", action="store_true", help="Skip the PNG picture.")
    args = ap.parse_args()

    cutoff = args.cutoff
    if cutoff is None:
        cutoff = 5.0 if args.method == "heavy" else 8.0

    chains_wanted = None
    if args.chains:
        chains_wanted = {c.strip() for c in args.chains.split(",") if c.strip()}

    prefix = args.prefix or os.path.splitext(os.path.basename(args.input))[0] or "rin"

    print(f"[1/4] loading structure: {args.input}")
    structure = load_structure(args.input, tmpdir=args.outdir)

    print(f"[2/4] selecting residues (model {args.model}, "
          f"chains {'all' if not chains_wanted else ','.join(sorted(chains_wanted))})")
    residues = collect_residues(structure, args.model, chains_wanted)
    if not residues:
        print("[error] no amino-acid residues found. Is this a protein "
              "structure? (nucleic-acid-only or all-hetero inputs yield no "
              "nodes.)", file=sys.stderr)
        sys.exit(2)

    print(f"[3/4] finding contacts: method={args.method}, cutoff={cutoff} A")
    if args.method == "heavy":
        labels, kept, coords, edges = build_by_heavy_atoms(
            residues, cutoff, args.min_seq_sep)
    else:
        labels, kept, coords, edges = build_by_representative(
            residues, args.method, cutoff, args.min_seq_sep)

    G = assemble_graph(labels, kept, coords, edges, args.weighted)

    print(f"[4/4] writing outputs to: {args.outdir}/")
    base = write_outputs(G, args.outdir, prefix, args.method, cutoff,
                         args.weighted, want_plot=not args.no_plot)

    print(f"\nDone. {G.number_of_nodes()} residues, {G.number_of_edges()} contacts.")
    print(f"Graph:     {base}.graphml")
    print(f"Adjacency: {base}_adjacency.npz  (+ {prefix}_nodes.csv for the order)")
    print(f"Summary:   {base}_summary.txt")


if __name__ == "__main__":
    main()
