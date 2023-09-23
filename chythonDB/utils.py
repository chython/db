# -*- coding: utf-8 -*-
#
#  Copyright 2023 Timur Gimadiev <timur.gimadiev@gmail.com>
#  Copyright 2023 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of chythonDB.
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
from datasketch import MinHashLSH, MinHashLSHEnsemble
from typing import Dict, Optional, Tuple
from datasketch.storage import _random_name

organic_set = {'B', 'C', 'N', 'O', 'Si', 'P', 'S', 'Se', 'F', 'Cl', 'Br', 'I'}


def validate_molecule(mol):
    if int(mol):
        return False
    for n, a in mol.atoms():
        if a.atomic_symbol not in organic_set:
            return False
    return True


class MinHashLSH(MinHashLSH):

    def _byteswap(self, hs):
        return int.from_bytes(hs.byteswap().data, 'big')

    def _hashed_byteswap(self, hs):
        return self.hashfunc(int.from_bytes(hs.byteswap().data, 'big'))


class MinHashLSHEnsemble(MinHashLSHEnsemble):
    class MinHashLSHEnsemble(MinHashLSHEnsemble):

        def __init__(
                self,
                threshold: float = 0.9,
                num_perm: int = 128,
                num_part: int = 16,
                m: int = 8,
                weights: Tuple[float, float] = (0.5, 0.5),
                storage_config: Optional[Dict] = None,
                prepickle: Optional[bool] = None,
        ) -> None:
            if threshold > 1.0 or threshold < 0.0:
                raise ValueError("threshold must be in [0.0, 1.0]")
            if num_perm < 2:
                raise ValueError("Too few permutation functions")
            if num_part < 1:
                raise ValueError("num_part must be at least 1")
            if m < 2 or m > num_perm:
                raise ValueError("m must be in the range of [2, num_perm]")
            if any(w < 0.0 or w > 1.0 for w in weights):
                raise ValueError("Weight must be in [0.0, 1.0]")
            if sum(weights) != 1.0:
                raise ValueError("Weights must sum to 1.0")
            self.threshold = threshold
            self.h = num_perm
            self.m = m
            rs = self._init_optimal_params(weights)
            # Initialize multiple LSH indexes for each partition
            storage_config = {"type": "dict"} if not storage_config else storage_config
            basename = storage_config.get("basename", _random_name(11))
            self.indexes = [
                dict(
                    (
                        r,
                        MinHashLSH(
                            num_perm=self.h,
                            params=(int(self.h / r), r),
                            storage_config=self._get_storage_config(
                                basename, storage_config, partition, r
                            ),
                            prepickle=prepickle,
                        ),
                    )
                    for r in rs
                )
                for partition in range(0, num_part)
            ]
            self.lowers = [None for _ in self.indexes]
            self.uppers = [None for _ in self.indexes]


__all__ = ['validate_molecule', 'MinHashLSH']
