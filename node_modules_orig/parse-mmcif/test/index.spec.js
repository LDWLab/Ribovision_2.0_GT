const { expect } = require('chai');
const { readFileSync } = require('fs');
const parseMmcif = require('../index');

describe('parsemmcif', () => {
  let mmcif;

  describe('when given 3AID', () => {
    beforeEach(() => {
      mmcif = readFileSync('./dat/3aid.cif', 'utf8');
    });

    it('reads the right number of atoms and reads the first atom correctly', () => {
      const { atoms } = parseMmcif(mmcif);
      const firstAtom = atoms[0];

      expect(atoms).to.have.lengthOf(1846);
      expect(firstAtom).to.have.property('id', 1);
      expect(firstAtom).to.have.property('type_symbol', 'N');
      expect(firstAtom).to.have.property('label_atom_id', 'N');
      expect(firstAtom).to.have.property('label_alt_id', '.');
      expect(firstAtom).to.have.property('label_comp_id', 'PRO');
      expect(firstAtom).to.have.property('label_asym_id', 'A');
      expect(firstAtom).to.have.property('label_entity_id', 1);
      expect(firstAtom).to.have.property('label_seq_id', 1);
      expect(firstAtom).to.have.property('pdbx_PDB_ins_code', '?');
      expect(firstAtom).to.have.property('Cartn_x', -2.555);
      expect(firstAtom).to.have.property('Cartn_y', 9.253);
      expect(firstAtom).to.have.property('Cartn_z', 34.411);
      expect(firstAtom).to.have.property('x', -2.555);
      expect(firstAtom).to.have.property('y', 9.253);
      expect(firstAtom).to.have.property('z', 34.411);
      expect(firstAtom).to.have.property('occupancy', 1.0);
      expect(firstAtom).to.have.property('B_iso_or_equiv', 30.6);
      expect(firstAtom).to.have.property('pdbx_formal_charge', '?');
      expect(firstAtom).to.have.property('auth_seq_id', 1);
      expect(firstAtom).to.have.property('auth_comp_id', 'PRO');
      expect(firstAtom).to.have.property('auth_asym_id', 'A');
      expect(firstAtom).to.have.property('auth_atom_id', 'N');
      expect(firstAtom).to.have.property('pdbx_PDB_model_num', 1);
    });
  });

  describe('when given super big file 6BO8', () => {
    beforeEach(() => {
      mmcif = readFileSync('./dat/6BO8.cif', 'utf8');
    });

    it('reads the right number of atoms and reads the first atom correctly', () => {
      console.log('justin do it');
      const { atoms } = parseMmcif(mmcif);
      const firstAtom = atoms[0];

      expect(atoms).to.have.lengthOf(19048);
      expect(firstAtom).to.have.property('id', 1);
      expect(firstAtom).to.have.property('type_symbol', 'N');
      expect(firstAtom).to.have.property('label_atom_id', 'N');
      expect(firstAtom).to.have.property('label_alt_id', '.');
      expect(firstAtom).to.have.property('label_comp_id', 'SER');
      expect(firstAtom).to.have.property('label_asym_id', 'A');
      expect(firstAtom).to.have.property('label_entity_id', 1);
      expect(firstAtom).to.have.property('label_seq_id', 28);
      expect(firstAtom).to.have.property('pdbx_PDB_ins_code', '?');
      expect(firstAtom).to.have.property('Cartn_x', 61.242);
      expect(firstAtom).to.have.property('Cartn_y', 81.517);
      expect(firstAtom).to.have.property('Cartn_z', 107.847);
      expect(firstAtom).to.have.property('x', 61.242);
      expect(firstAtom).to.have.property('y', 81.517);
      expect(firstAtom).to.have.property('z', 107.847);
      expect(firstAtom).to.have.property('occupancy', 1.0);
      expect(firstAtom).to.have.property('B_iso_or_equiv', 183.72);
      expect(firstAtom).to.have.property('pdbx_formal_charge', '?');
      expect(firstAtom).to.have.property('auth_seq_id', 28);
      expect(firstAtom).to.have.property('auth_comp_id', 'SER');
      expect(firstAtom).to.have.property('auth_asym_id', 'A');
      expect(firstAtom).to.have.property('auth_atom_id', 'N');
      expect(firstAtom).to.have.property('pdbx_PDB_model_num', 1);
    });
  });
});
