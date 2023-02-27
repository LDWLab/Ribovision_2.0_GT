# Advanced features

## Synchronization of navigation between the panels
Navigation between all panels is synchronized. Hovering over an alignment position highlights the corresponding residue in the Topology and MolStar viewers. Reversely, hovering over a residue in the Topology or MolStar viewers highlights the alignment and the other structural viewer. When the amino-acid frequencies are shown, hovering over the datapoints on the graph highlights the current species in the alignment viewer.

## Guide
ProteoVision offers an interactive guide that demonstrates the steps a user can take to fully utilize ProteoVision capabilities. The guide always starts on the users first visit and can be launched at any time with the **Help** button. The user can step through the guide using the keyboard arrow keys or the **Next** and **Previous** buttons on the guide pop-ups. Ending the guide with **Skip tour** button will erase the current session and reset the viewports.

## User upload mode
ProteoVision supports the upload of custom alignments. The user can upload any fasta format alignment through the **User upload** menu and calculate amino-acid frequencies and mapping data from it. The user can select a PDB ID to visualize a structure and map the calculated data from their alignment.

## Selecting a structure for visualization and mapping in user upload mode
After uploading an alignment, ProteoVision will use the first sequence of the alignment to perform a BLAST search of the [PDB database](https://www.ebi.ac.uk/Tools/common/tools/help). BLAST results are filtered by E-value lower than 10<sup>-5</sup> and are used to populate a dropdown menu in the PDB input field. The dropdown menu is searchable, and shows filtered results depending on the user input. After a BLAST is complete the polymers associated with the polymer selection box will be filtered by the BLAST results.
A BLAST search can take a long time, so ProteoVision displays a message that BLAST is running under the PDB input field (“BLASTing available PDBs”). While waiting the user can input any 4 letter PDB ID in the PDB input field the fetched polymers for that PDB ID will not be filtered. When the BLAST search is complete the message will change to indicate completion (“Completed BLAST for similar PDBs”).

## DESIRE-API
ProteoVision uses an API service which provides data about rProtein nomenclature, sequences, alignments, and annotations. Furthermore, it provides data about species phylogeny. The API is available [here]( https://proteovision.chemistry.gatech.edu/desire-api/).
