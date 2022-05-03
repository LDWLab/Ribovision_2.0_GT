# Basic Navigation 

## Overview
ProteoVision is operated via the main Navigation panel for selection/filtering and contains three main applets:
- Alignment Viewer for representation of multiple sequence alignments; 
- PDB topology Viewer for depiction of protein secondary structures; 
- MolStar viewer for visualization of three-dimensional structures. 

Additionally, non-mappable statistical data (e.g. amino acids frequencies) are visualized in a separate window using Plotly applet. Detailed description of all functions for ProteoVision is listed in this documentation. ProteoVision also contains an optional interactive guided tour with a brief description of each functional element.

## Selecting phylogenetic group(s) 
To retrieve a specified alignment from the DESIRE database, follow the steps below: 

1. Click on the dropdown menu **Select a phylogenetic group**; the three dropdown nodes are labeled with the three domains of life: Eukarya, Bacteria, and Archaea.
2. Click on the left-side drop-down arrow of any phylogenetic group in the dropdown menu; the phylogenetic subgroups of that group become visible. For example, upon clicking the left-side arrow of “Bacteria,” the following groups are made visible to the user: “Bacteroidetes,” “Chlorobi,” “Cyanobacteria,” “Proteobacteria,” etc. This step is recursive.
3. Add a phylogenetic group to the list of phylogenetic groups for which the alignment will be retrieved. This list is displayed at the top of the drop-down structure. Note that when zero groups are selected this field displays the text **Select a phylogenetic group**. Multiple groups can be selected.
4. The user can also directly search for a phylogenetic group by typing in its name.

## Selecting a protein alignment
Upon selection of a set of phylogenetic groups (See Selecting phylogenetic group(s)), a list of compatible protein alignments is pulled from the DESIRE database in a new dropdown menu. The list of alignments is updated after every new addition in the phylogenetic browser, the alignments list is filtered to show alignments that contain all selected phylogenetic groups. For example, if a user selects Archaea, Bacteria and Eukarya, only universal rProtein alignments will be shown since only universal proteins are present in all three TOL branches. The user selects a protein alignment from the available ones and the alignment is displayed in the alignment viewer; the rows of the alignment are phylogenetic groups, and the columns of the alignment are arranged according to individual amino acid indices. 
Some rProteins can have more than one alignment. For example, aL18 is a protein present in Archaea and Eukarya, and there are three possible alignments for it: one that includes only archaeal sequences (aL18), one that includes only Eukaryotic ones (eL18) and one that includes both (aeL18). The alignment will still be filtered to include only species selected in the phylogenetic browser. For example, if the user selects Eukarya from the phylogenetic browser and aeL18 from the alignment menu, they will fetch an alignment that was built from archaeal and eukaryotic sequences but is truncated down to show only eukaryotic sequences.

## Selecting a structure for visualization and mapping
Once an alignment of a specified protein is selected, the option to select its structure for visualization and mapping becomes available. Selecting a structure is a two-step process:

1. First the user selects a 3D structure available from PDBe. Structures from the PDBe are filtered by the polymer of the selected alignment. In DESIRE mode the filtering is done by connecting to [riboXYZ](https://ribosome.xyz) and retrieving PDB IDs with matching rProtein nomenclature. In user upload mode the filtering is done by using the first sequence of the provided alignment to BLAST the RCSB database and retrieve highly scoring PDB IDs. In both modes the user can search the filtered PDB IDs by typing the desired PDB ID or the species name. In both modes the user is free to type in any PDB ID with length 4 and ProteoVision will try to find relevant chains from the input. For user convenience, the first three PDB ID options are never filtered and are those for the cytosolic ribosomes of E. coli (4V9D), P. furiosus (4V6U), and H. sapiens (4UG0). It is possible to encounter Eukaryotic ribosomal structures even if the user has selected a bacterial protein; these structures will be non-cytoplasmic (e.g., plastid or mitochondrial).
2. Once a PDB ID is selected, the available chains of that structure are filtered by the polymer identity present in the selected alignment. No filtering is done when using a custom alignment and all chains of the structure are available. The user can select one chain which will load the topology and MolStar viewers. 

## Creating a structure-alignment mapping
To ensure the match between the selected MSA and the sequence from a selected structure and to properly visualize the data, we compute an alignment between the sequences of the structure and the alignment. This is done internally on the server using the mafft program with the –addfull option (Katoh; [10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010) and it does not require any action from the user. However, if the structural sequence has extra positions compared to the alignment, ProteoVision displays an error message located between the alignment and 2D/3D viewers. The message warns the user how many positions failed to map properly.

## Selecting attribute data to map
The user has an option to map either calculated mapping data (see Advanced features) or custom data (supplied by a user) onto the 2D and 3D viewers. Available data attributes may be selected from a dropdown menu in the lower right corner of the topology viewer. There is also an option in the sidebar to upload custom mapping data. This data should be uploaded as a .csv; an example file is provided on the ProteoVision site. Once the custom data has been uploaded, it will be added as a mapping option in the dropdown menu of the PDB topology viewer applet. 
