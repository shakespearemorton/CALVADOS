"""Read single chains from PDB and PDBx/mmCIF files by author chain ID."""

import warnings

import numpy as np
from MDAnalysis import AtomGroup, Universe
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.lib.util import convert_aa_code
from numpy.typing import NDArray
from openmm import app
from openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader

from .inputmodels import InputPath

RNA_RESNAMES = {"A", "C", "G", "U"}
DNA_RESNAMES = {"DA", "DC", "DG", "DT", "DU", "DI"}


def _cif_author_chain_map(cif: str) -> dict[str, str]:
    """Map the chain IDs OpenMM assigned to the file's author chain IDs.

    OpenMM's PDBxFile switches to ``label_asym_id`` when it has more unique
    values than ``auth_asym_id`` (e.g. ligands and waters sharing the protein's
    author chain). Each label chain belongs to exactly one author chain.
    """
    with open(cif) as f:
        data: list = []
        PdbxReader(f).read(data)
    atoms = data[0].getObj("atom_site")
    auth = atoms.getAttributeIndex("auth_asym_id")
    label = atoms.getAttributeIndex("label_asym_id")
    if auth == -1 or label == -1:
        return {}
    rows = atoms.getRowList()
    if len({row[label] for row in rows}) > len({row[auth] for row in rows}):
        return {row[label]: row[auth] for row in rows}
    return {}


def load_universe(path: InputPath) -> Universe:
    """Load a PDB or CIF structure whose ``chainIDs`` are author chain IDs."""
    path = str(path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if path.lower().endswith(".cif"):
            u = Universe(app.pdbxfile.PDBxFile(path))
            chain_map = _cif_author_chain_map(path)
            if chain_map:
                u.atoms.chainIDs = [chain_map.get(c, c) for c in u.atoms.chainIDs]
        else:
            u = Universe(path)
    return u


def residue_kind(resname: str) -> str | None:
    """Classify a residue name as protein, rna, dna, or non-polymer (None)."""
    resname = resname.strip()
    if resname in RNA_RESNAMES:
        return "rna"
    if resname in DNA_RESNAMES:
        return "dna"
    if len(resname) == 3:
        try:
            convert_aa_code(resname)
        except ValueError:
            return None
        return "protein"
    return None


def polymer_chains(u: Universe) -> dict[str, str]:
    """Return the protein and RNA chains of a structure in file order.

    DNA chains and non-polymer residues (waters, ions, ligands) are skipped
    with a warning.
    """
    kinds: dict[str, set[str]] = {}
    skipped: dict[str, set[str]] = {}
    for res in u.residues:
        chain = str(res.atoms[0].chainID)
        kind = residue_kind(res.resname)
        if kind is None:
            skipped.setdefault(chain, set()).add(res.resname)
        else:
            kinds.setdefault(chain, set()).add(kind)

    chains = {}
    for chain, chain_kinds in kinds.items():
        if chain_kinds == {"dna"}:
            warnings.warn(f"Skipping DNA chain {chain!r}: DNA is not supported.")
        elif len(chain_kinds) > 1:
            raise ValueError(f"Chain {chain!r} mixes residue types {sorted(chain_kinds)}.")
        else:
            chains[chain] = chain_kinds.pop()
    for chain, resnames in skipped.items():
        warnings.warn(f"Skipping non-polymer residues {sorted(resnames)} of chain {chain!r}.")
    return chains


def select_chain(u: Universe, chain: str) -> AtomGroup:
    """Return the atoms of the protein/RNA residues of one author chain."""
    indices = [
        idx for idx, res in enumerate(u.residues)
        if str(res.atoms[0].chainID) == chain
        and residue_kind(res.resname) in ("protein", "rna")
    ]
    return u.residues[indices].atoms


def protein_sequence(ag: AtomGroup) -> str:
    """Return the one-letter sequence of the protein residues in ``ag``."""
    return "".join(convert_aa_code(res.resname) for res in ag.residues)


def positions(
    ag: AtomGroup, com_mask: NDArray[np.bool_], origin: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Return one position (nm) per residue relative to ``origin`` (nm).

    Residues where ``com_mask`` is True use their centre of mass, the rest
    their CA atom.
    """
    pos = ag.select_atoms("name CA").positions
    if com_mask.any():
        coms = np.array([res.atoms.center_of_mass() for res in ag.residues])
        pos = np.where(com_mask[:, None], coms, pos)
    return pos / 10.0 - origin

def go_pairs_from_structure(
    xyz_i: NDArray[np.float64],
    xyz_j: NDArray[np.float64],
    cutoff: float = 0.95,
) -> list[tuple[int, int, float]]:
    """Return Go contacts between the beads of two molecules.

    ``xyz_i`` and ``xyz_j`` (nm) must share a frame, e.g. two chains of a
    complex. Contacts are bead pairs closer than ``cutoff`` (nm), returned as
    0-based bead indices with their native distance (nm).
    """
    dist = distance_array(xyz_i, xyz_j)
    return [
        (int(i), int(j), float(dist[i, j]))
        for i, j in zip(*np.nonzero(dist < cutoff))
    ]
