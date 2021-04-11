# Saving
## Saving the alignment (fasta & image)
The alignment retrieved from the DESIRE database can be downloaded in fasta format with the **Download alignment** button. 
The current viewport of the alignment can also be saved as a png image with the **Download alignment image** dropdown menu. 
The user can select whether to download the currently visible region of the alignment or the entire alignment.


## Saving a CD-HIT report
When using custom alignment ProteoVision runs [CD-HIT](http://weizhongli-lab.org/cd-hit/) to ensure there is no overrepresentation in the sequences. The user can download the complete CD-HIT report from the CD-HIT dropdown menu above the alignment.
## Saving secondary structure image (svg) 
The secondary structure image may be downloaded as a .svg file by clicking the “**S**” button in the bottom right corner of the topology viewer. 

## Saving 3D structure image (png) 
The image of the 3D structure may be downloaded as a .png file or viewed in browser by clicking the **Screenshot/State** button, which appears as a wheel in the upper right corner of the 3D viewer. Upon clicking this button, the user may choose to either keep the default white background or make the background transparent by turning transparency to **On**. They may also choose to either include or exclude 3D axes which show the orientation of the protein.

## Saving computed data as csv file or pymol script 
Calculated properties for a given alignment can be saved as .csv file or as a pymol script from the **Download mapped data** button. For the .csv file the first column is labeled Index and indicates the residue number, the rest of the columns hold the calculated attributes for each alignment column mapped on the current structure indices. The pymol script fetches the selected PDB ID, extracts the selected chain, creates separate objects for each of the calculated properties, and colors them the same way they are colored in the 3D viewer.

## Saving frequencies (png) 
The amino acid frequency plot can be saved as a .png file by hovering over the image and clicking on the **camera icon** in the top center. 

## Saving a ProteoVision session
At any point, the user can save their progress with the **Save session** button. This will download a .json file that holds information about the currently loaded alignment, structure, and custom mapping data. The session file does not save information about the masking and truncation ranges.
