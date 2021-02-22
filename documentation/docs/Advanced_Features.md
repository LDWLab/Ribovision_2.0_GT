# Advanced features
## Frequencies 
Once an alignment has been selected, the user can click the **show amino-acid frequencies** button to display amino acid frequencies for that alignment. 
Amino acid frequencies are calculated by processing the alignment fasta using Biopython, and the resulting data is displayed as a faceted boxplot through the Plotly graphing library.
Every point on an amino acid frequency plot represents the frequency of a single amino acid for a single species within the alignment. The species and amino acid frequency associated with each point can be viewed by hovering over the point.
If a user chooses from the **Select secondary structure dropdown**, the amino acid frequencies will be calculated for the specified secondary structure within the polypeptide. 
Users can choose to display frequencies from helix residues, coil residues, or strand residues.

