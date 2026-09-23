#!/usr/bin/env python3
"""rin_builder.py - protein structure -> residue interaction network (a graph).

Input: a PDB/mmCIF file or a 4-character PDB ID (downloaded). Output: the graph
as GraphML + sparse adjacency (.npz) + node/edge tables + a summary + a picture.
Nodes are amino-acid residues; edges join residues in contact. Handles the messy
real cases: multi-model NMR, multiple chains, waters/hetero, alternate
conformations, missing atoms, insertion codes, and modified residues (e.g. MSE
-> MET). Pure Python, cross-platform.

    python rin_builder.py 1A8O.pdb            # or a PDB ID
    python rin_builder.py file.cif --method heavy --chains A
"""
import argparse, os, sys, warnings
import numpy as np, networkx as nx
from scipy import sparse
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.simplefilter("ignore", PDBConstructionWarning)

# modified residues mapped to their standard parent so they aren't dropped
MODIFIED_TO_STANDARD = {
    "MSE": "MET", "SEP": "SER", "TPO": "THR", "PTR": "TYR", "HYP": "PRO",
    "PCA": "GLN", "CSO": "CYS", "MLY": "LYS", "M3L": "LYS", "KCX": "LYS",
    "CME": "CYS", "CSD": "CYS", "OCS": "CYS", "FME": "MET",
}


def load_structure(inp, tmpdir="."):
    """Parse a structure from a file path or fetch it by 4-character PDB ID."""
    if os.path.isfile(inp):
        parser = (MMCIFParser(QUIET=True)
                  if os.path.splitext(inp)[1].lower() in (".cif", ".mmcif")
                  else PDBParser(QUIET=True))
        return parser.get_structure("structure", inp)
    token = inp.strip()
    if len(token) == 4 and token.isalnum():
        from Bio.PDB import PDBList
        pdbl = PDBList()
        path = pdbl.retrieve_pdb_file(token, pdir=tmpdir, file_format="pdb")
        if not (path and os.path.isfile(path)):
            path = pdbl.retrieve_pdb_file(token, pdir=tmpdir, file_format="mmCif")
            return MMCIFParser(QUIET=True).get_structure(token, path)
        return PDBParser(QUIET=True).get_structure(token, path)
    raise FileNotFoundError(f"'{inp}' is neither a readable file nor a 4-char PDB ID.")


def standard_resname(residue):
    """Standard 3-letter name if this residue is (or maps to) an amino acid, else None."""
    name = residue.get_resname().strip().upper()
    if name in MODIFIED_TO_STANDARD:
        return MODIFIED_TO_STANDARD[name]
    return name if is_aa(name, standard=True) else None


def pick_representative_atom(residue, method):
    """(name, coord) for the residue's representative atom, or None. 'cb' uses
    CB, falling back to CA (glycine has no CB); 'ca' uses CA."""
    def coord(a):
        return residue[a].get_coord() if a in residue else None
    if method == "ca":
        c = coord("CA")
        return ("CA", c) if c is not None else None
    c = coord("CB")
    if c is not None:
        return ("CB", c)
    c = coord("CA")
    return ("CA", c) if c is not None else None


def collect_residues(structure, model_index, chains_wanted):
    """Ordered list of accepted amino-acid residues in the chosen model."""
    models = list(structure.get_models())
    if not models:
        raise ValueError("structure contains no models.")
    if model_index >= len(models):
        raise ValueError(f"requested model {model_index} but only {len(models)} present.")
    out = []
    for chain in models[model_index]:
        cid = chain.get_id().strip() or "_"
        if chains_wanted and cid not in chains_wanted:
            continue
        for res in chain:
            std = standard_resname(res)
            if std is None:
                continue  # water, ligand, ion, nucleic acid, etc.
            _, resseq, icode = res.get_id()
            out.append({"chain": cid, "resseq": resseq, "icode": (icode or "").strip(),
                        "resname": std, "orig_resname": res.get_resname().strip(),
                        "residue": res})
    return out


def node_label(info):
    """Stable node id like 'A:57' or 'A:57B' (insertion code)."""
    return f"{info['chain']}:{info['resseq']}{info['icode']}"


def _seq_far_enough(a, b, min_seq_sep):
    """True unless a and b are same-chain and closer than min_seq_sep in sequence."""
    if min_seq_sep <= 0 or a["chain"] != b["chain"]:
        return True
    return abs(a["resseq"] - b["resseq"]) >= min_seq_sep


def build_by_representative(residues, method, cutoff, min_seq_sep):
    """Contacts via one representative atom per residue (CB/CA) within cutoff."""
    labels, coords, kept = [], [], []
    for info in residues:
        rep = pick_representative_atom(info["residue"], method)
        if rep and rep[1] is not None:
            labels.append(node_label(info))
            coords.append(rep[1])
            kept.append(info)
    coords = np.asarray(coords, dtype=float)
    edges = []
    if len(coords) >= 2:
        for i, j in cKDTree(coords).query_pairs(r=cutoff):
            if _seq_far_enough(kept[i], kept[j], min_seq_sep):
                edges.append((labels[i], labels[j], float(np.linalg.norm(coords[i] - coords[j]))))
    return labels, kept, coords, edges


def build_by_heavy_atoms(residues, cutoff, min_seq_sep):
    """Contacts via the minimum distance between any heavy atoms of two residues."""
    labels = [node_label(info) for info in residues]
    atoms, owner = [], []
    for idx, info in enumerate(residues):
        for atom in info["residue"]:
            if atom.get_name().startswith("H") or (atom.element or "").strip().upper() == "H":
                continue
            atoms.append(atom)
            owner.append(idx)
    edges = {}
    if len(atoms) >= 2:
        atom_owner = {id(a): owner[k] for k, a in enumerate(atoms)}
        for a, b in NeighborSearch(atoms).search_all(cutoff, level="A"):
            ia, ib = atom_owner[id(a)], atom_owner[id(b)]
            if ia == ib or not _seq_far_enough(residues[ia], residues[ib], min_seq_sep):
                continue
            key = (min(ia, ib), max(ia, ib))
            d = float(np.linalg.norm(a.get_coord() - b.get_coord()))
            if key not in edges or d < edges[key]:
                edges[key] = d
    coords = []
    for info in residues:
        rep = pick_representative_atom(info["residue"], "ca")
        coords.append(rep[1] if rep and rep[1] is not None else (np.nan,) * 3)
    coords = np.asarray(coords, dtype=float)
    return labels, residues, coords, [(labels[i], labels[j], d) for (i, j), d in edges.items()]


def assemble_graph(labels, kept, coords, edges, weighted):
    """NetworkX graph with node attributes and edge distances/weights."""
    G = nx.Graph()
    for k, info in enumerate(kept):
        xyz = (0.0, 0.0, 0.0)
        if coords is not None and not np.any(np.isnan(coords[k])):
            xyz = tuple(float(v) for v in coords[k])
        G.add_node(labels[k], chain=info["chain"], resseq=int(info["resseq"]),
                   icode=info["icode"], resname=info["resname"],
                   orig_resname=info["orig_resname"], x=xyz[0], y=xyz[1], z=xyz[2])
    for u, v, d in edges:
        G.add_edge(u, v, distance=round(d, 3),
                   weight=round((1.0 / d) if (weighted and d > 0) else 1.0, 6))
    return G


def write_outputs(G, outdir, prefix, method, cutoff, weighted, want_plot):
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, prefix)
    nx.write_graphml(G, base + ".graphml")
    nodes = list(G.nodes(data=True))

    with open(base + "_nodes.csv", "w") as fh:
        fh.write("index,node,chain,resseq,icode,resname,x,y,z\n")
        for i, (n, at) in enumerate(nodes):
            fh.write(f"{i},{n},{at['chain']},{at['resseq']},{at['icode']},"
                     f"{at['resname']},{at['x']:.3f},{at['y']:.3f},{at['z']:.3f}\n")
    with open(base + "_edges.csv", "w") as fh:
        fh.write("source,target,distance,weight\n")
        for u, v, at in G.edges(data=True):
            fh.write(f"{u},{v},{at['distance']},{at['weight']}\n")

    order = [n for n, _ in nodes]
    A = nx.to_scipy_sparse_array(G, nodelist=order,
                                 weight=("weight" if weighted else None), format="csr")
    sparse.save_npz(base + "_adjacency.npz", A)

    degs = [d for _, d in G.degree()]
    with open(base + "_summary.txt", "w") as fh:
        fh.write("Residue interaction network summary\n-----------------------------------\n")
        fh.write(f"contact method     : {method}\n")
        fh.write(f"distance cutoff (A) : {cutoff}\n")
        fh.write(f"edge weighting      : {'1/distance' if weighted else 'binary'}\n")
        fh.write(f"residues (nodes)    : {G.number_of_nodes()}\n")
        fh.write(f"contacts (edges)    : {G.number_of_edges()}\n")
        fh.write(f"graph density       : {nx.density(G) if G.number_of_nodes() > 1 else 0:.4f}\n")
        fh.write(f"connected comps     : {nx.number_connected_components(G) if G.number_of_nodes() else 0}\n")
        if degs:
            fh.write(f"mean/min/max degree : {np.mean(degs):.2f} / {min(degs)} / {max(degs)}\n")
        fh.write(f"chains included     : {', '.join(sorted({at['chain'] for _, at in nodes}))}\n")

    if want_plot and G.number_of_nodes():
        try:
            import matplotlib; matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            pos = nx.spring_layout(G, seed=1)
            deg = np.array(degs, dtype=float)
            plt.figure(figsize=(9, 9))
            nx.draw_networkx_edges(G, pos, alpha=0.2, width=0.5)
            nx.draw_networkx_nodes(G, pos, node_size=25 + 4 * deg, node_color=deg, cmap="viridis")
            plt.axis("off"); plt.title(f"{prefix}: {G.number_of_nodes()} residues, {G.number_of_edges()} contacts")
            plt.tight_layout(); plt.savefig(base + "_network.png", dpi=140); plt.close()
        except Exception as e:
            print(f"[warn] plot skipped: {e}", file=sys.stderr)
    return base


def main():
    ap = argparse.ArgumentParser(description="Protein structure -> residue interaction network.")
    ap.add_argument("input", help="PDB/mmCIF file path, or a 4-character PDB ID.")
    ap.add_argument("--method", choices=["cb", "ca", "heavy"], default="cb",
                    help="cb: CB-CB (CA for Gly). ca: CA-CA. heavy: min heavy-atom distance.")
    ap.add_argument("--cutoff", type=float, default=None, help="Cutoff (A). Default 8 cb/ca, 5 heavy.")
    ap.add_argument("--chains", default=None, help="Chains to include (default all).")
    ap.add_argument("--model", type=int, default=0, help="Model index for NMR files.")
    ap.add_argument("--min-seq-sep", type=int, default=0,
                    help="Drop same-chain contacts closer than this in sequence.")
    ap.add_argument("--weighted", action="store_true", help="1/distance edge weights.")
    ap.add_argument("--outdir", default="rin_output")
    ap.add_argument("--prefix", default=None)
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()

    cutoff = a.cutoff if a.cutoff is not None else (5.0 if a.method == "heavy" else 8.0)
    chains = {c.strip() for c in a.chains.split(",") if c.strip()} if a.chains else None
    prefix = a.prefix or os.path.splitext(os.path.basename(a.input))[0] or "rin"

    print(f"[1/4] loading structure: {a.input}")
    structure = load_structure(a.input, tmpdir=a.outdir)
    print(f"[2/4] selecting residues (model {a.model}, chains {'all' if not chains else ','.join(sorted(chains))})")
    residues = collect_residues(structure, a.model, chains)
    if not residues:
        sys.exit("[error] no amino-acid residues found (nucleic-acid-only or all-hetero input).")
    print(f"[3/4] finding contacts: method={a.method}, cutoff={cutoff} A")
    parts = (build_by_heavy_atoms(residues, cutoff, a.min_seq_sep) if a.method == "heavy"
             else build_by_representative(residues, a.method, cutoff, a.min_seq_sep))
    G = assemble_graph(*parts, a.weighted)
    print(f"[4/4] writing outputs to: {a.outdir}/")
    base = write_outputs(G, a.outdir, prefix, a.method, cutoff, a.weighted, not a.no_plot)
    print(f"\nDone. {G.number_of_nodes()} residues, {G.number_of_edges()} contacts.")
    print(f"Graph: {base}.graphml  |  Adjacency: {base}_adjacency.npz  |  Summary: {base}_summary.txt")


if __name__ == "__main__":
    main()
