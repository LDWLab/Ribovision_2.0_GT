# Advanced features
## Frequencies 
Once an alignment has been selected, the user can click the **show amino-acid frequencies** button to display amino acid frequencies for that alignment. 
Amino acid frequencies are calculated by processing the alignment fasta using Biopython, and the resulting data is displayed as a faceted boxplot through the Plotly graphing library.
Every point on an amino acid frequency plot represents the frequency of a single amino acid for a single species within the alignment. The species and amino acid frequency associated with each point can be viewed by hovering over the point.
If a user chooses from the **Select secondary structure dropdown**, the amino acid frequencies will be calculated for the specified secondary structure within the polypeptide. 
Users can choose to display frequencies from helix residues, coil residues, or strand residues.

## Masking 
Once a polymer has been selected, the user may select to **mask/unmask 2D and 3D residues**. This allows for coloration of only selected residue ranges specified by the user. All other residues will be colored in white and will have hovering functions disabled. The overall structure of the protein will still be visible.

## Range selection 
The user also has the option to **cut/uncut 2D and 3D residues**. This allows the user to view only the portion of the protein specified by the entered residues. The rest of the protein structure will be removed rather than colored white as in the masking feature. 

## Synchronization of navigation between the panels
Navigation between all panels is synchronized. Hovering over an alignment position highlights the corresponding residue in the Topology and MolStar viewers. Reversely, hovering over a residue in the Topology or MolStar viewers highlights the alignment and the other structural viewer. When the amino-acid frequencies are shown, hovering over the datapoints on the graph highlights the current species in the alignment viewer.

## Guide
ProteoVision offers an interactive guide that demonstrates the steps a user can take to fully utilize ProteoVision capabilities. The guide always starts on the users first visit and can be launched at any time with the **Help** button. The user can step through the guide using the keyboard arrow keys or the **Next** and **Previous** buttons on the guide pop-ups. Ending the guide with **Skip tour** button will erase the current session and reset the viewports.

## User upload mode

### Custom alignments
ProteoVision supports the upload of custom alignments. The user can upload any fasta format alignment through the **User upload** menu and calculate amino-acid frequencies and mapping data from it. The user can select a PDB ID or upload a custom PDB file to visualize a structure and map the calculated data from their alignment.

#### CD-HIT analysis
When using custom alignment ProteoVision runs [CD-HIT](http://weizhongli-lab.org/cd-hit/) to ensure there is no overrepresentation in the sequences. By default ProteoVision uses 90% identity clustering threshold. Sequences with greater similarity are automatically removed from the uploaded alignment. The user can elect to use their original alignment by selecting the option from the CD-HIT dropdown menu above the alignment.

### Custom structures
ProteoVision supports the upload of custom PDB structures. The web server supports PDB files with only a single chain. The user can upload PDB format structure through the **Upload a custom PDB** button (available after uploading a custom alignment). ProteoVision connects the structural sequence to the alignment, displays the tertiary structure in the Mol* viewer, generates a topology representation with [Pro-Origami](http://munk.cis.unimelb.edu.au/pro-origami/), and displays the topology with the Topology Viewer.

## Selecting a structure for visualization and mapping in user upload mode
After uploading an alignment, ProteoVision will use the first sequence of the alignment to perform a BLAST search of the [PDB database](https://www.ebi.ac.uk/Tools/common/tools/help). BLAST results are filtered by E-value lower than 10<sup>-5</sup> and coverage of the query sequence greater than 75%. The BLAST results are used to populate a dropdown menu in the PDB input field. The dropdown menu is searchable, and shows filtered results depending on the user input. After a BLAST is complete the polymers associated with the polymer selection box will be filtered by the BLAST results.
A BLAST search can take a long time, so ProteoVision displays a message that BLAST is running under the PDB input field (“BLASTing available PDBs”). While waiting the user can input any 4 letter PDB ID in the PDB input field the fetched polymers for that PDB ID will not be filtered. When the BLAST search is complete the message will change to indicate completion (“Completed BLAST for similar PDBs”).

## DESIRE-API
ProteoVision uses an API service which provides data about rProtein nomenclature, sequences, alignments, and annotations. Furthermore, it provides data about species phylogeny. The API is available [here]( https://proteovision.chemistry.gatech.edu/desire-api/).
