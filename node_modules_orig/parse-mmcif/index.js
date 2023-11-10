const ATOM_NAME = 'ATOM';

const FIELDS = [
  'type',
  'id',
  'type_symbol',
  'label_atom_id',
  'label_alt_id',
  'label_comp_id',
  'label_asym_id',
  'label_entity_id',
  'label_seq_id',
  'pdbx_PDB_ins_code',
  'Cartn_x',
  'Cartn_y',
  'Cartn_z',
  'occupancy',
  'B_iso_or_equiv',
  'pdbx_formal_charge',
  'auth_seq_id',
  'auth_comp_id',
  'auth_asym_id',
  'auth_atom_id',
  'pdbx_PDB_model_num',
];

const INTEGER_FIELDS = new Set([
  'id',
  'label_entity_id',
  'label_seq_id',
  'auth_seq_id',
  'pdbx_PDB_model_num',
]);

const FLOAT_FIELDS = new Set([
  'Cartn_x',
  'Cartn_y',
  'Cartn_z',
  'occupancy',
  'B_iso_or_equiv',
]);

/**
 * Parses the given mmCIF string into json
 * @param {String} mmcif
 * @returns {Object}
 */
module.exports = function parseMmcif(mmcif) {
  const mmcifLines = mmcif.split('\n');
  const atoms = [];

  // Iterate each line looking for atoms
  mmcifLines.forEach((mmcifLine) => {
    if (mmcifLine.substr(0, 4) === ATOM_NAME) {
      // http://mmcif.wwpdb.org/docs/tutorials/content/atomic-description.html
      let fieldIndex = 0;
      const atom = mmcifLine.split(' ').reduce((atom, value) => {
        if (value === '') {
          return atom;
        }

        // Properly parse the value
        const field = FIELDS[fieldIndex];
        if (INTEGER_FIELDS.has(field)) {
          atom[field] = parseInt(value);
        } else if (FLOAT_FIELDS.has(field)) {
          atom[field] = parseFloat(value);
        } else {
          atom[field] = value;
        }

        // Alias x, y, z coords
        if (field === 'Cartn_x') {
          atom.x = atom[field];
        } else if (field === 'Cartn_y') {
          atom.y = atom[field];
        } else if (field === 'Cartn_z') {
          atom.z = atom[field];
        }

        fieldIndex++;
        return atom;
      }, {});
      atoms.push(atom);
    }
  });

  return {
    atoms,
  };
}

