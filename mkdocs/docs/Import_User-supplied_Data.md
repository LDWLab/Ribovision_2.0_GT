# Import User-supplied Data
## Import an external dataset for a pre-selected alignment
Once an alignment and a structure have been selected, the option becomes available to upload custom data for mapping onto the selected structure for visualization in the topology and MolStar viewers. 
The data should be supplied in a .csv (comma-separated values) file format. The first row of the csv file contains the headers for each column. An Index column should always be indicated. 
All other additional columns require a unique header definition. The Index column has the residue number to which the user-supplied data will be mapped. The rest of the columns have the data that will be mapped in the form of numerical values. 
Once a correct csv file has been uploaded, the selected structure will be colored in the topology and MolStar viewers according to the values in the csv file. 
Different columns in the csv file will appear as different **Annotations** in the dropdown menu in the lower right corner of the topology viewer. 

## Upload an external multiple sequence alignment. 
The **User Upload** mode allows the user to use all the ProteoVision features with an external MSA. An alignment in a .fasta file format must be selected and then uploaded with the **Upload alignment** button. Once an alignment is uploaded, it will be displayed in the Alignment viewer. The steps for structure selection, mapping, attribute calculation, and saving are the same as previously described.

## Upload an external PDB structure file. 
The **User Upload** mode allows the user to use all the ProteoVision features with an external structure. A structural model in a .pdb file format containing a single chain must be selected and then uploaded with the **Upload a custom PDB** button (available after uploading a custom alignment). Once a structure is uploaded, ProteoVision will process it and display 2D and 3D representations. The steps for mapping, attribute calculation, and saving are the same as previously described.

## Import a saved ProteoVision session file. 
At any time, the user can upload a previously saved ProteoVision session file with the **Load session** button. This will restore the progress at time of saving for selected/uploaded alignment, calculated/uploaded data attributes, and selected structures in the 2D and 3D viewers. Currently the sessions do not recover masking or truncation ranges.
