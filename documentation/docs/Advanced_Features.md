# Advanced features

## Synchronization of navigation between the panels
Navigation between all panels is synchronized. Hovering over an alignment position highlights the corresponding residue in the Topology and MolStar viewers. Reversely, hovering over a residue in the Topology or MolStar viewers highlights the alignment and the other structural viewer. When the amino-acid frequencies are shown, hovering over the datapoints on the graph highlights the current species in the alignment viewer.

## Show protein contacts

Ribovision 2.0 offers visualization of rProtein contacts with a selected rRNA molecule. Processing is performed by NeighborSearch module of Biopython on the fly upon selection an desired RNA chain in the given PDB complex. The list of proteins in contact with a given RNA molecule (the Protein Contact List) apprears in the lower part of the Main Naviagation panel after the processing is fininshed.  All contacts within a cut off of 3.5 A are reported. Thus, LSU may contain a sub-set of rProteins from the SSU that are intersubunit bridges, and vice versa. The protein contacts are visualized in 2D and 3D viewer applets upon clicking on desired sub-set of protein in the Protein Contact List. In addition, the 3D structures of the selected rProteins appear in the 3D applet.  


## Show modified nucleotides

Ribovision 2.0 offers visualization of modified nuclotides in a selected rRNA molecule. If a selected RNA chanin contain modifiedd nucleotides, the abbreviated modification types (e.g. 1MG, PSU, 5MU etc) will appear in the List of Modified Nuclotides beneath of Protein Contact List in the Main Navigation Pannel. The modified residues will be highlighted in the 2D and 3D applets  upon selecting the desired modification type from the List of Modified Nuclotides. All nucleototides of a given modification type will be highlighted in the same color. Selection of multiple modification types will higlight modified nucleotides in distinct colors (per modification type).
 


## DESIRE-API
ProteoVision uses an API service which provides data about rProtein nomenclature, sequences, alignments, and annotations. Furthermore, it provides data about species phylogeny. The API is available [here]( https://ribovision2.chemistry.gatech.edu/desire-api/).
