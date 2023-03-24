# Advanced features

## Synchronization of navigation between the panels
Navigation between all panels is synchronized. Hovering over an alignment position highlights the corresponding residue in the Topology and MolStar viewers. Reversely, hovering over a residue in the Topology or MolStar viewers highlights the alignment and the other structural viewer. When the amino-acid frequencies are shown, hovering over the datapoints on the graph highlights the current species in the alignment viewer.

## Show protein contacts

Ribovision 2.0 offers visualization of rProtein contacts with a selected rRNA molecule. Processing is performed by NeighborSearch module of Biopython on the fly upon selection an desired RNA chain in the given PDB complex. The list of proteins in contact with a given RNA molecule (the Protein Contact List) apprears in the lower part of the Main Naviagation panel after the processing is fininshed.  All contacts within a cut off of 3.5 A are reported. Thus, LSU may contain a sub-set of rProteins from the SSU that are intersubunit bridges, and vice versa. The protein contacts are visualized in 2D and 3D viewer applets upon clicking on desired sub-set of protein in the Protein Contact List. In addition, the 3D structures of the selected rProteins appear in the 3D applet.  


## Show modified nucleotides

Ribovision 2.0 offers visualization of modified nuclotides in a selected rRNA molecule. If a selected RNA chanin contain modifiedd nucleotides, the abbreviated modification types (e.g. 1MG, PSU, 5MU etc) will appear in the List of Modified Nuclotides beneath of Protein Contact List in the Main Navigation Pannel. The modified residues will be highlighted in the 2D and 3D applets  upon selecting the desired modification type from the List of Modified Nuclotides. All nucleototides of a given modification type will be highlighted in the same color. Selection of multiple modification types will higlight modified nucleotides in distinct colors (per modification type).


## User upload data
RiboVision 2.0  supports the visualization of user supplied data. Once an alignment and a structure have been selected, the option to upload custom data for mapping onto the selected structure for visualization in the topology and MolStar viewers becomes available in the Main Navigation Panel. The data should be supplied in a .csv (comma-separated values) file format. The first row of the csv file contains the headers for each column. An Index column should always be indicated. All other additional columns require a unique header definition. The Index column has the residue number to which the user-supplied data will be mapped. The rest of the columns have the data that will be mapped in the form of numerical values. Once a correct csv file has been uploaded, the selected structure will be colored in the topology and MolStar viewers according to the values in the csv file. Different columns in the csv file will appear as different Annotations in the Select Data dropdown menu in the lower right corner of the RNA topology viewer.

## Selecting a structure for visualization and mapping in user upload mode
After uploading an alignment, ProteoVision will use the first sequence of the alignment to perform a BLAST search of the [PDB database](https://www.ebi.ac.uk/Tools/common/tools/help). BLAST results are filtered by E-value lower than 10<sup>-5</sup> and are used to populate a dropdown menu in the PDB input field. The dropdown menu is searchable, and shows filtered results depending on the user input. After a BLAST is complete the polymers associated with the polymer selection box will be filtered by the BLAST results.
A BLAST search can take a long time, so ProteoVision displays a message that BLAST is running under the PDB input field (“BLASTing available PDBs”). While waiting the user can input any 4 letter PDB ID in the PDB input field the fetched polymers for that PDB ID will not be filtered. When the BLAST search is complete the message will change to indicate completion (“Completed BLAST for similar PDBs”).



## DESIRE-API
ProteoVision uses an API service which provides data about rProtein nomenclature, sequences, alignments, and annotations. Furthermore, it provides data about species phylogeny. The API is available [here]( https://proteovision.chemistry.gatech.edu/desire-api/).
