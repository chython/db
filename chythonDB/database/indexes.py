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
import json
from types import MethodType
from datasketch import MinHash
from datasketch.lsh import _optimal_param
from pony.orm import Required, PrimaryKey, IntArray, db_session, select
from tqdm import tqdm

from . import config, db
from .config import Entity, Config
from ..utils import MinHashLSH, MinHashLSHEnsemble
from uuid import uuid4
from collections import namedtuple

RequestPack = namedtuple('RequestPack', ['request', 'postprocess', 'prefetch_map'])


class MoleculeSimilarityIndex(Entity):
    band = Required(int)
    key = Required(int, size=64)
    records = Required(IntArray)
    PrimaryKey(band, key)


class CGRSimilarityIndex(Entity):
    band = Required(int)
    key = Required(int, size=64)
    records = Required(IntArray)
    PrimaryKey(band, key)


class MoleculeContainmentIndex(Entity):
    partition = Required(int)
    model = Required(int)
    band = Required(int)
    key = Required(int, size=64)
    records = Required(IntArray)
    PrimaryKey(partition, model, band, key)

def fingerprintidx_create(self):
    print("creation of index for fingerprints")
    self.execute("DROP INDEX if exists idx_moleculestructure;")
    self.execute("CREATE INDEX idx_moleculestructure ON MoleculeStructure USING GIST(fingerprint gist__intbig_ops);")
    self.execute("DROP INDEX if exists idx_cgrfingerprint;")
    self.execute("CREATE INDEX idx_cgrfingerprint ON cgr USING GIST(fingerprint gist__intbig_ops);")
    print("creation of index for fingerprint length")
    self.execute("DROP INDEX if exists idx_cgr_fingerprint_len;")
    self.execute("CREATE INDEX idx_cgr_fingerprint_len ON cgr(fingerprint_len);")
    self.execute("DROP INDEX if exists idx_fingerprint_len;")
    self.execute("CREATE INDEX idx_fingerprint_len ON moleculestructure(fingerprint_len);")


def molecule_similatiry_idx_create(self):
    num_permute = config.lsh_num_permute or 64
    threshold = config.lsh_threshold or 0.7
    lsh = MinHashLSH(threshold=threshold, num_perm=num_permute, hashfunc=hash)
    b, r = _optimal_param(threshold, num_permute, 0.5, 0.5)
    hashranges = [(i * r, (i + 1) * r) for i in range(b)]
    with db_session:
        self.execute(f"""DELETE FROM MoleculeSimilarityIndex WHERE 1=1""")
        self.execute(f"""DELETE FROM Config c WHERE c.key='hashranges' """)
        self.commit()
        Config(key="hashranges", value=json.dumps(hashranges))
        print("Creation of LSH")
        for idx, fingerprint in tqdm(self.execute(f'SELECT id, fingerprint FROM MoleculeStructure')):
            h = MinHash(num_perm=num_permute, hashfunc=hash)
            h.update_batch(fingerprint)
            lsh.insert(idx, h, check_duplication=False)
        print("Uploading LSH tables into DB")
        for n, ht in tqdm(enumerate(lsh.hashtables, 1)):
            for key, value in ht._dict.items():
                self.insert(MoleculeSimilarityIndex, band=n, key=key, records=list(value))
            self.commit()

def molecule_containment_idx_create(self):
    num_permute = config.lsh_num_permute or 64
    threshold = config.lsh_threshold or 0.7
    lsh = MinHashLSHEnsemble(threshold=threshold, num_perm=num_permute)
    with db_session:
        self.execute(f"""DELETE FROM MoleculeContainmentIndex WHERE 1=1""")
        self.execute(f"""DELETE FROM Config c WHERE c.key='hashranges' """)
        self.commit()
        print("Creation of LSH")
        tmp = []
        for idx, fingerprint in tqdm(self.execute(f'SELECT id, fingerprint FROM MoleculeStructure')):
            h = MinHash(num_perm=num_permute, hashfunc=hash)
            h.update_batch(fingerprint)
            tmp.append((idx, h, len(fingerprint)))
        lsh.index(tmp)
        print("Uploading LSH tables into DB")
        hashranges = {}
        for idx, model in lsh.indexes[0].items():
            hashranges[idx] = model.hashranges
        Config(key="lsh_ensemble_hashranges", value=json.dumps(hashranges))
        for partition, i in enumerate(lsh.indexes):
            for lsh_mod, (y, z) in enumerate(i.items()):
                for band, g in enumerate(z.hashtables):
                    for key, records in g._dict.items():
                        self.insert(MoleculeContainmentIndex(partition=partition, model=lsh_mod, band=band, key=key,
                                                             records=records))


def cgr_similatiryidx_create(self):
    num_permute = config.cgr_lsh_num_permute or 64
    threshold = config.cgr_lsh_threshold or 0.7
    lsh = MinHashLSH(threshold=threshold, num_perm=num_permute, hashfunc=hash)
    b, r = _optimal_param(threshold, num_permute, 0.5, 0.5)
    hashranges = [(i * r, (i + 1) * r) for i in range(b)]
    with db_session:
        self.execute(f"""DELETE FROM CGRSimilarityIndex WHERE 1=1""")
        self.execute(f"""DELETE FROM Config c WHERE c.key='cgr_hashranges' """)
        self.commit()
        Config(key="cgr_hashranges", value=json.dumps(hashranges))
        print("Creation of LSH for CGR")
        for idx, fingerprint in tqdm(self.execute(f'SELECT id, fingerprint FROM cgr')):
            h = MinHash(num_perm=num_permute, hashfunc=hash)
            h.update_batch(fingerprint)
            lsh.insert(idx, h, check_duplication=False)
        print("Uploading CGR LSH tables into DB")
        for n, ht in tqdm(enumerate(lsh.hashtables, 1)):
            for key, value in ht._dict.items():
                self.insert(CGRSimilarityIndex, band=n, key=key, records=list(value))
            self.commit()


class CursorHolder:
    __slots__ = ('function', 'cursor', 'buffer', 'prefetch_class', 'prefetch_position', 'prefetch_fields')

    def __init__(self, request_pack: RequestPack):
        self.function = request_pack.postprocess
        self.prefetch_class, self.prefetch_position, self.prefetch_fields = request_pack.prefetch_map
        self.cursor = db.get_connection().cursor(str(uuid4()))# withhold=True)
        self.cursor.execute(request_pack.request)
        self.buffer = self.loader()

    def loader(self):
        if res := self.cursor.fetchone():
            if self.prefetch_fields:
                i = res[self.prefetch_position]
                select(x for x in self.prefetch_class if x.id == i).prefetch(*self.prefetch_fields)[:]
            yield res
        else:
            return

        c = 10
        while True:
            if res := self.cursor.fetchmany(c):
                if self.prefetch_fields:
                    ids_to_prefetch = [x[self.prefetch_position] for x in res]
                    select(x for x in self.prefetch_class if x.id in ids_to_prefetch).prefetch(*self.prefetch_fields)[:]
                yield from res
                if c < 1000:
                    c *= 10
            else:
                return

    def __iter__(self):
        return self

    def __next__(self):
        while True:
            if result := self.function(next(self.buffer)):
                if result == "stop":
                    raise StopIteration("Similarity limit reached")
                return result

    #def __del__(self):
    #    if not self.cursor.closed:
    #        self.cursor.close()

def unbind(self):
    self.provider = self.schema = None

db.create_sim_index = MethodType(molecule_similatiry_idx_create, db)
db.create_containment_index = MethodType(molecule_containment_idx_create, db)
db.create_fing_index = MethodType(fingerprintidx_create, db)
db.unbind = MethodType(unbind, db)
db.create_cgr_sim_index = MethodType(cgr_similatiryidx_create, db)


__all__ = ['MoleculeSimilarityIndex', 'molecule_similatiry_idx_create', 'molecule_containment_idx_create',
           'RequestPack', 'CursorHolder']
