# Advanced features

## Data in the 2D RNA topology viewer
Once the secondary structure is fetched from the PDBe API (or generated by R2DT) and mapped to the MSA, users can i) display various representations of the 2D structure (nucleotides, contour lines, and circles); visualize various types of base-pairs according to Leontis-Westhof notation (all or nested only); iii) map precomputed properties (e.g. Shannon entropy or TwinCons) onto the secondary structure (the selected data will be synchronously mapped onto the 3D representation as well).
Users also have ability to save the resulting secondary structure representations (along with the mapped data) as SVG files.


## Synchronization of navigation between the panels
Navigation between all panels is synchronized. Hovering over an alignment position highlights the corresponding residue in the topology and Mol* viewers. Correspondingly, hovering over a residue in the Topology or Mol* viewers highlights the alignment and the other structural viewer. 


## Show protein contacts
Ribovision 2.0 offers visualization of rProtein contacts with a selected rRNA molecule. Processing is performed by the NeighborSearch module of Biopython on the fly upon selection of a desired RNA chain in the given PDB complex. The list of proteins in contact with a given RNA molecule (the Protein Contact List) apprears in the lower part of the Main Naviagation panel after the processing is finished.  All contacts within a cut off of 3.5 A are reported. Thus, LSU contacts may contain a subset of rProteins from the SSU that are intersubunit bridges, and vice versa. The protein contacts are visualized in the 2D and 3D viewer applets upon clicking on a desired subset of proteins in the Protein Contact List. In addition, the 3D structures of the selected rProteins appear in the 3D applet.  A tool tip will display additional information (protein name and its chain in the given PDB) upon hovering the mouse over a highlighted nucleotide in the 2D applet. 


## Show modified nucleotides
Ribovision 2.0 offers visualization of modified nuclotides in a selected rRNA molecule. If a selected RNA chain contains modified nucleotides, the abbreviated modification types (e.g. 1MG, PSU, 5MU etc) will appear in the List of Modified Nuclotides beneath the Protein Contact List in the Main Navigation Panel. The modified residues will be highlighted in the 2D and 3D applets upon selecting the desired modification type from the List of Modified Nuclotides. All nucleototides of a given modification type will be highlighted in the same color. Selection of multiple modification types will higlight modified nucleotides in distinct colors (per modification type). A tool tip will display additional information (modification name) upon hovering the mouse on a highlighted nucleotide in the 2D applet. 


## User upload data
RiboVision 2.0  supports the visualization of user supplied data. Once an alignment and a structure have been selected, the option to upload custom data for mapping onto the selected structure for visualization in the topology and Mol* viewers becomes available in the Main Navigation Panel. The data should be supplied in a .csv (comma-separated values) file format. The first row of the csv file contains the headers for each column. An Index column should always be indicated. All other additional columns require a unique header definition. The Index column has the residue number to which the user-supplied data will be mapped. The rest of the columns have the data that will be mapped in the form of numerical values. Once a correct csv file has been uploaded, the selected structure will be colored in the topology and Mol* viewers according to the values in the csv file. Different columns in the csv file will appear as different Annotations in the Select Data dropdown menu in the lower right corner of the RNA topology viewer.


## DESIRE-API
ProteoVision uses an API service which provides data about rProtein nomenclature, sequences, alignments, and annotations. Furthermore, it provides data about species phylogeny. The API is available [here]( https://ribovision2.chemistry.gatech.edu/desire-api/).
