# Basic Navigation 

## Overview
ProteoVision is operated via the main Navigation panel for selection/filtering and contains three main applets:
- Alignment Viewer for representation of multiple sequence alignments; 
- RNA topology Viewer for depiction of RNA secondary structures layouts; 
- MolStar viewer for visualization of three-dimensional structures. 


## Selecting phylogenetic group(s) 
To retrieve a specified alignment from the DESIRE database, follow the steps below: 

1. Click on the dropdown menu **Select a phylogenetic group**; the three dropdown nodes are labeled with the three domains of life: Eukarya, Bacteria, and Archaea.
2. Click on the left-side drop-down arrow of any phylogenetic group in the dropdown menu; the phylogenetic subgroups of that group become visible. For example, upon clicking the left-side arrow of “Bacteria,” the following groups are made visible to the user: “Bacteroidetes,” “Chlorobi,” “Cyanobacteria,” “Proteobacteria,” etc. This step is recursive.
3. Add a phylogenetic group to the list of phylogenetic groups for which the alignment will be retrieved. This list is displayed at the top of the drop-down structure. Note that when zero groups are selected this field displays the text **Select a phylogenetic group**. Multiple groups can be selected.
4. The user can also directly search for a phylogenetic group by typing in its name.

## Selecting an RNA alignment
Upon selection of a set of phylogenetic groups (See Selecting phylogenetic group(s)), a list of molecules (LSU or SSU) appears, and once a desired group is selected, compatible RNA alignments is pulled from the RiboVision2.0 database in a new dropdown menu. The list of alignments is updated after every new addition in the phylogenetic browser, the alignments list is filtered to show alignments that contain all selected phylogenetic groups. For example, if a user selects Archaea, Bacteria and Eukarya, only universal RNA alignments will be shown. The user selects an RNA alignment from the available ones and the alignment is displayed in the alignment viewer; the rows of the alignment are phylogenetic groups, and the columns of the alignment are arranged according to individual nucleotide indices. .

## Selecting a structure for visualization and mapping
Once an alignment of a specified protein is selected, the option to select its structure for visualization and mapping becomes available. Selecting a structure is a two-step process:

1. First the user selects a 3D structure available from PDBe. Structures from the PDBe are filtered by the polymer of the selected alignment. Filtered PDB IDs are available in a dropdown menu when using the DESIRE database. When using a custom alignment, no filtering is done on PDBs and the user can write in any 4 letter PDB ID. For user convenience, the first three PDB ID options are never filtered and are those for the cytosolic ribosomes of E. coli (7K00), P. furiosus (4V6U), and H. sapiens (4UG0). 
2. Once a PDB ID is selected, the available chains of that structure are filtered by the polymer identity present in the selected alignment. No filtering is done when using a custom alignment and all chains of the structure are available. The user can select one chain which will load the topology and MolStar viewers. 

## Creating a structure-alignment mapping
To ensure the match between the selected MSA and the sequence from a selected structure and to properly visualize the data, we compute an alignment between the sequences of the structure and the alignment. This is done internally on the server using the mafft program with the –addfull option (Katoh; [10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)) and it does not require any action from the user. However, if the structural sequence has extra positions compared to the alignment, RiboVision 2.0 displays an error message located between the alignment and 2D/3D viewers. The message warns the user how many positions failed to map properly.

## Selecting attribute data to map
The user has an option to map either calculated mapping data (see Advanced features) or custom data (supplied by a user) onto the 2D and 3D viewers. Available data attributes may be selected from a dropdown menu in the lower right corner of the topology viewer. There is also an option in the sidebar to upload custom mapping data. This data should be uploaded as a .csv; an example file is provided on the RiboVision2.0 site. Once the custom data has been uploaded, it will be added as a mapping option in the dropdown menu of the PDB topology viewer applet. 
