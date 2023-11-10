# parse-mmcif
A small utility for parsing MMCIF files into useable JSON.

## Getting Started

    npm install --save parse-pdb

```
const parseMmcif = require('parse-mmcif');
const { readFileSync } = require('fs');

const pdbString = readFileSync('./3aid.cif', 'utf8');

const parsed = parsePdb(pdbString);

console.log(parsed.atoms);
/*
[ {
    type: 'ATOM',
    id: 1,
    type_symbol: 'N',
    label_atom_id: 'N',
    label_alt_id: '.',
    label_comp_id: 'PRO',
    label_asym_id: 'A',
    label_entity_id: 1,
    label_seq_id: 1,
    pdbx_PDB_ins_code: '?',
    Cartn_x: -2.555,
    x: -2.555,
    Cartn_y: 9.253,
    y: 9.253,
    Cartn_z: 34.411,
    z: 34.411,
    occupancy: 1,
    B_iso_or_equiv: 30.6,
    pdbx_formal_charge: '?',
    auth_seq_id: 1,
    auth_comp_id: 'PRO',
    auth_asym_id: 'A',
    auth_atom_id: 'N',
    pdbx_PDB_model_num: 1,
  },
  ...1845 others
]
*/
```

## JSON Format
The output json is an object containing arrays of each structure keyed on record name, according to the [pdb spec](http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html).

```
  atoms:
    type: string
    id: integer
    type_symbol: string
    label_atom_id: string
    label_alt_id: string
    label_comp_id: string
    label_asym_id: string
    label_entity_id: integer
    label_seq_id: integer
    pdbx_PDB_ins_code: string
    Cartn_x: float
    x: -float
    Cartn_y: float
    y: float
    Cartn_z: float
    z: float
    occupancy: integer
    B_iso_or_equiv: float
    pdbx_formal_charge: string
    auth_seq_id: integer
    auth_comp_id: string
    auth_asym_id: string
    auth_atom_id: string
    pdbx_PDB_model_num: integer
```

## Upcoming
I plan to add support for residues, chains, and anything else someone might need.  Pull requests welcome!

## License
MIT.  See LICENSE file.
