# Import User-supplied Data
## Import an external dataset for a pre-selected alignment
Once an alignment and a structure have been selected, the option becomes available to upload custom data for mapping onto the selected structure for visualization in the topology and MolStar viewers. 
The data should be supplied in a .csv (comma-separated values) file format. The first row of the csv file contains the headers for each column. An Index column should always be indicated. 
All other additional columns require a unique header definition. The Index column has the residue number to which the user-supplied data will be mapped. The rest of the columns have the data that will be mapped in the form of numerical values. 
Once a correct csv file has been uploaded, the selected structure will be colored in the topology and MolStar viewers according to the values in the csv file. 
Different columns in the csv file will appear as different **Annotations** in the dropdown menu in the lower right corner of the topology viewer. 
