# User-upload Mode
The **User Upload** mode allows the user to use most of the features of RiboVision 2.0 with an external MSA and 3D structure.

## Upload an external multiple sequence alignment. 
A custom alignment (or a single RNA sequence) in a .fasta file format must be chosen and uploaded with the **Upload alignment** button. Once an alignment is uploaded, it will be displayed in the Alignment viewer.

## Upload an external 3D structure.
The **User Upload** mode requires the user to upload the the structure file  in one of the two formats (CIF or PDB). 
1) The CIF cormat may contain multiple chains but is requred to contain a certain fields described below:
   
a) following ids:
-auth_seq_id, auth_comp_id, label_entity_id, auth_asym_id;

b) block containing: 

-loop_

_entity_poly.entity_id
_entity_poly.type
_entity_poly.nstd_linkage
_entity_poly.nstd_monomer
_entity_poly.pdbx_seq_one_letter_code
_entity_poly.pdbx_seq_one_letter_code_can
_entity_poly.pdbx_strand_id
_entity_poly.pdbx_target_identifier

This block is required by RiboVision 2.0 code to extract the complete RNA sequence from the cif file.

c)  block containing:
-loop_ 
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
_pdbx_poly_seq_scheme.ndb_seq_num
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.auth_seq_num
_pdbx_poly_seq_scheme.pdb_mon_id
_pdbx_poly_seq_scheme.auth_mon_id
_pdbx_poly_seq_scheme.pdb_strand_id
_pdbx_poly_seq_scheme.pdb_ins_code
_pdbx_poly_seq_scheme.hetero

This block is required to map the full (genomic) RNA sequence onto nucleotides resolved in the structure (even if the all nucleotides are resolved in the 3D structure). RiboVision 2.0 code takes an advantage of the provided cif dictionary and performs the mapping automatically.

In addition to uploading of the CIF file, the users must specify an entity ID of the desired RNA chain. This will be done in a field form beneath "Upload CIF file" button.

2) The PDB file is restricted to a single RNA chain.
In addition to uploading the 3D structure in the PDB file, the users must supply a full genetic sequence of the corresponding RNA chain. This requirement ensures that the secondary structure is generated from the complete sequence (and not from the resolved fragments).    Upon uploading the genomic RNA sequence, it will be appended to the MSA (as full sequence ) by Mafft. This sequence will be excluded from all statistics (e.g. nucleotide frequency).

The hallmark of RiboVision 2.0 is generation of the  2D RNA structure on the fly using a template based R2DT algorithm. This process may take several minutes (depending on the size of the supplied RNA chain). The available RNA templates include but not limited to the LSU and SSU rRNAs, tRNAs, RNAseP, microRNAs. Please refer to R2DT documentation for additional details and current capabilities of [R2DT](https://docs.r2dt.bio/en/latest/). 

Upon completion of R2DT job, the  secondary structure will appear in the RNA topology viewer, and will be interactively linked to the MSA  and the uploaded 3D structure. Base pairs are currently derived from the R2DT layout. No 3D derived base pairing option are currenly suported. Thus, only canonical (cWW and Wobble)  base pairings are generated in the User-upload mode.


## Import an external dataset for a pre-selected alignment
Once an alignment and a structure have been selected, the option becomes available to upload custom data for mapping onto the selected structure for visualization in the topology and MolStar viewers. 
The data should be supplied in a .csv (comma-separated values) file format. The first row of the csv file contains the headers for each column. An Index column should always be indicated. 
All other additional columns require a unique header definition. The Index column has the residue number to which the user-supplied data will be mapped. The rest of the columns have the data that will be mapped in the form of numerical values. 
Once a correct csv file has been uploaded, the selected structure will be colored in the topology and MolStar viewers according to the values in the csv file. 
Different columns in the csv file will appear as different **Custom Data** in the **Select Data** dropdown menu  of the RNA topology viewer. 




