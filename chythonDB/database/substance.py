# -*- coding: utf-8 -*-
#
#  Copyright 2021 Timur Gimadiev <timur.gimadiev@gmail.com>
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of CGRdb.
#
#  CGRdb is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from .config import Entity
from . import molecule
from ..utils import validate_molecule

from collections.abc import Iterable as abcIterable
from functools import cached_property, reduce
from operator import or_
from typing import Tuple, Iterable, Union, Dict, List, Optional
from chython import MoleculeContainer
from pony.orm import PrimaryKey, Required, Set, Optional as PonyOptional, Json


class Substance(Entity):
    id = PrimaryKey(int, auto=True)
    reactions = Set('ReactionSubstance')
    components = Set('SubstanceStructure')
    name = PonyOptional(str)

    def __init__(self, substance: Optional[Iterable[Tuple[MoleculeContainer, Optional[float]]]] = None, /, name=None):
        """
        :param substance:
        """
        if name:
            super().__init__(name=name)
        else:
            super().__init__()
        if substance is not None:
            if not isinstance(substance, abcIterable):
                raise ValueError("Iterable type consisted of tuples should be provided")
            substance = list(substance)
            for mol, frac in substance:
                if not isinstance(mol, MoleculeContainer):
                    raise ValueError('First element of tuple should be CGRtools.MoleculeContainer')
                if frac is not None and not isinstance(frac, float):
                    raise TypeError('molar_fraction should be of float type')

            for mol, frac in substance:
                sao = mol.smiles_atoms_order
                if validate_molecule(mol):
                    if (ms := molecule.MoleculeStructure.get(signature=bytes(mol))) is None:
                        structure = molecule.Molecule(mol)
                        mapping = None
                    else:
                        structure = ms.molecule
                        mapping = dict(zip(ms._fast_mapping, sao))
                elif (structure := molecule.NonOrganic.get(signature=bytes(mol))) is None:
                    structure = molecule.NonOrganic(mol)
                    mapping = None
                else:
                    mapping = dict(zip(structure._fast_mapping, sao))
                SubstanceStructure(structure, self, molar_fraction=frac, mapping=mapping)

    @cached_property
    def structure(self) -> MoleculeContainer:
        return reduce(or_, (x.structure for x in self.components))


class SubstanceStructure(Entity):
    id = PrimaryKey(int, auto=True)
    molar_fraction = PonyOptional(float)
    _mapping = PonyOptional(Json, column='mapping', nullable=True)
    substance = Required('Substance')
    molecule = PonyOptional('Molecule')
    non_organic = PonyOptional('NonOrganic')

    def __init__(self, structure: Union['molecule.Molecule', 'molecule.NonOrganic'], substance: Substance, /, *,
                 molar_fraction: Optional[float] = None,
                 mapping: Union[Dict[int, int], List[Tuple[int, int]], None] = None):
        if isinstance(mapping, dict):
            if not all(isinstance(x, int) for x in mapping.items() for x in x):
                raise TypeError('Mapping is dict that contains following structure: '
                                '{int (atom of substance) : int (atom of NonOrganic or Molecule)}')
            mapping = list(mapping.items())
        elif isinstance(mapping, (tuple, list)):
            if not all(isinstance(x, (tuple, list)) and len(x) == 2 and all(isinstance(x, int) for x in x)
                       for x in mapping):
                raise TypeError('For the tuple pairs of integers following constructions should be applied: '
                                '[(int (atom of substance), int (atom of NonOrganic or Molecule)),...]')
        elif mapping is not None:
            raise TypeError('Mapping of Substance to NonOrganic or Molecule expected')

        if isinstance(structure, molecule.Molecule):
            mol = structure
            non_organic = None
        elif isinstance(structure, molecule.NonOrganic):
            mol = None
            non_organic = structure
        else:
            raise TypeError('At least one either Molecule or NonOrganic required')
        super().__init__(molar_fraction=molar_fraction, molecule=mol, _mapping=mapping,
                         non_organic=non_organic, substance=substance)

    @cached_property
    def structure(self):
        if self.molecule:
            if self.mapping:
                return self.molecule.canonic_structure.structure.remap(self.mapping, copy=True)
            return self.molecule.canonic_structure.structure
        else:
            if self.mapping:
                return self.non_organic.structure.remap(self.mapping, copy=True)
            return self.non_organic.structure

    @cached_property
    def mapping(self):
        return dict(self._mapping) if self._mapping else {}


__all__ = ['SubstanceStructure', 'Substance']
