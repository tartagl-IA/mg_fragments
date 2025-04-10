"""Filters for molecules."""

import json
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from rdkit.Chem import Mol, MolFromSmarts, rdMolDescriptors

import settings


@dataclass
class ReactivePattern:
    """Class to represent a reactive pattern."""

    name: str
    smarts: str
    mol: Optional[Mol] = None

    def to_dict(self) -> dict:
        """Convert the ReactivePattern to a JSON serializable dictionary."""
        return {
            "name": self.name,
            "smarts": self.smarts,
        }


def get_reactive_pattern_list() -> list[ReactivePattern]:
    """Load reactive patterns from a JSON file.

    Returns:
        list[ReactivePattern]: List of reactive patterns.
    """
    with open(settings.REACTIVE_PATTERN_FILE, "r") as f:
        reactive_pattern_list = [ReactivePattern(**p) for p in json.load(f)]
    for reactive_pattern in reactive_pattern_list:
        reactive_pattern.mol = MolFromSmarts(reactive_pattern.smarts)
    return reactive_pattern_list


class Flexibility(str, Enum):
    """Enum to represent the flexibility of a molecule.

    Attributes:
        RIGID: Rigid molecule.
        FLEXIBLE: Flexible molecule.
    """

    RIGID = "rigid"
    FLEXIBLE = "flexible"

    @classmethod
    def values(cls) -> list[str]:
        """Get all possible values of the Flexibility enum.

        Returns:
            list[str]: List of all possible values.
        """
        return [e.value for e in cls]


def mol_reactive(
    mol: Mol, reactive_pattern_list: None | list[ReactivePattern] = None
) -> bool:
    """Check if a molecule contains electrophilic groups using precompiled SMARTS.

    Args:
        mol (Mol): Molecule to check.

    Returns:
        bool: True if the molecule contains reactive groups, False otherwise.
    """
    if reactive_pattern_list is None:
        reactive_pattern_list = get_reactive_pattern_list()
    return any(mol.HasSubstructMatch(pattern.mol) for pattern in reactive_pattern_list)


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
    mol: Mol, flexibility: Flexibility, max_rotable_bond_flexible: None | int = None
) -> bool:
    """Check if the molecule is flexible or rigid.

    Args:
        mol (Mol): RDKit molecule object.
        flexibility (Flexibility): Flexibility type.
        max_rotable_bond_flexible (int): Maximum number of rotatable bonds for flexible molecules.
            Default None (unset).

    Returns:
        bool: True if the molecule is flexible or rigid according to the specified criteria.
    """
    rb_num = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if not max_rotable_bond_flexible:
        max_rotable_bond_flexible = 1
    if (flexibility == Flexibility.RIGID and rb_num == 0) or (
        flexibility == Flexibility.FLEXIBLE and rb_num == max_rotable_bond_flexible
    ):
        return True
    return False
