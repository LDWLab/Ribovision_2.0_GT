# User-upload Mode
The **User Upload** mode allows the user to use most of the features of RiboVision 2.0 with an external MSA and 3D structure.

## Upload an external multiple sequence alignment. 
An alignment in a .fasta file format must be selected and then uploaded with the **Upload alignment** button. Once an alignment is uploaded, it will be displayed in the Alignment viewer.

## Upload an external 3D structure.
The **User Upload** mode allows the user to upload the the structure file  in the .PDB format. Currently, the pdb file must cotain a single RNA chain. Upon uploading the RNA sequence from the PDB file in extracted and appended to the supplied MSA by Mafft. Additionally, the 2D RNA structure is generated on the fly using a template based R2DT algorithm. This process may take several minutes. Upon completion of R2DT job, the  secondary structure will appear in the RNA topology viewer, and will be interactively linked to the MSA  and the uploaded 3D structure. Base pairs are currently derived from the generated R2DT layout. No 3D derived base pairing option are currenly suported. Thus, only canonical (cWW and Wobble)  base pairings are generated in the User-upload mode. The nucleotides unresolved in the provided 3D PDB file are ignored in the current iimplementation of the algorithm. The supprt for the CIF files that contain information  about the full polymer sequences of RNA chains is currently in progress.


## Import an external dataset for a pre-selected alignment
Once an alignment and a structure have been selected, the option becomes available to upload custom data for mapping onto the selected structure for visualization in the topology and MolStar viewers. 
The data should be supplied in a .csv (comma-separated values) file format. The first row of the csv file contains the headers for each column. An Index column should always be indicated. 
All other additional columns require a unique header definition. The Index column has the residue number to which the user-supplied data will be mapped. The rest of the columns have the data that will be mapped in the form of numerical values. 
Once a correct csv file has been uploaded, the selected structure will be colored in the topology and MolStar viewers according to the values in the csv file. 
Different columns in the csv file will appear as different **Annotations** in the dropdown menu in the lower right corner of the topology viewer. 





