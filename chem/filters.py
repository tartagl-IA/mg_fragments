"""Filters for molecules."""
import json
from enum import Enum
from rdkit.Chem import Mol, MolFromSmarts, rdMolDescriptors

# REACTIVE_PATTERN_SMARTS_DICT = [
#     {
#         "name": "Michael acceptor",
#         "smarts": "C=O[C]=[C]",
#     },
#     {
#         "name": "Aldehydes",
#         "smarts": "[CX3H1](=O)[#6]",
#     },
#     {
#         "name": "Acyl halides",
#         "smarts": "[CX3](=[OX1])[F,Cl,Br,I]",
#     },
#     {
#         "name": "Isocyanates",
#         "smarts": "[CX3](=O)N=C=O",
#     },
#     {
#         "name": "Sulfonyl halides",
#         "smarts": "[SX4](=O)(=O)[Cl,Br,I]",
#     },
#     {
#         "name": "Epoxides",
#         "smarts": "[C;R][O;R][C;R]",
#     },
#     {
#         "name": "Aziridines",
#         "smarts": "[#6]-1-[#6]-[#7]-1",
#     },
#     {
#         "name": "Thiiranes",
#         "smarts": "[C;R][S;R][C;R]",
#     }
# ]

def get_reactive_patterns():
    with open("chem/reactive_patterns.json", "r") as f:
        patterns = json.load(f)
    print(patterns)
    for pattern in patterns:
        pattern["mol"] = MolFromSmarts(pattern["smarts"])

    return patterns

REACTIVE_PATTERN_SMARTS_DICT_CACHE = [
    pattern["mol"] for pattern in get_reactive_patterns()
]

class Flexibility(str, Enum):
    RIGID = "rigid"
    FLEXIBLE = "flexible"
    
    @classmethod
    def values(cls) -> list[str]:
        return [e.value for e in cls]


def mol_reactive(mol: Mol, mol_patterns = None) -> bool:
    """Check if a molecule contains electrophilic groups using precompiled SMARTS.

    Args:
        mol (Mol): Molecule to check.

    Returns:
        bool: True if the molecule contains reactive groups, False otherwise.
    """
    if not mol_patterns:
        patterns = REACTIVE_PATTERN_SMARTS_DICT_CACHE
    return any(
        mol.HasSubstructMatch(pattern) for pattern in mol_patterns
    )

def mol_dimension_range(mol: Mol, min_atoms: int = 0, max_atoms: int = 0) -> bool:
    """Check if the dimension of the molecule is between the acceptance values.

    Args:
        mol (Mol): RDKit molecule object.
        min_atoms (int): minimum number of atoms. Default 0.
        max_atoms (int): maximum number of atoms. Default 0 (unset).

    Returns:
        bool: True if the molecule is between the defined range. False otherwise.
    """
    if max_atoms == 0:
        return min_atoms <= mol.GetNumAtoms()
    return min_atoms <= mol.GetNumAtoms() <= max_atoms

def mol_flexibility(
    mol: Mol, 
    flexibility: Flexibility,
    max_rotable_bond_flexible: None | int = None
) -> bool:
    rb_num = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if not max_rotable_bond_flexible:
        max_rotable_bond_flexible = 1
    if (
        flexibility == Flexibility.RIGID and rb_num == 0 or
        flexibility == Flexibility.FLEXIBLE and rb_num == max_rotable_bond_flexible
    ):
        return True
    return False
