# RiboVision 2.0 Data
## Phylogeny (DESIRE)
The subset of 152 species from the [DESIRE](https://doi.org/10.1093/molbev/msy101) (Sparse and Efficient Representation of Extant Biology) 
database was organized into a phylogenetic browser using a tree topology from the [Banfiled lab](https://doi.org/10.1038/nmicrobiol.2016.48).

## RNA Sequences
Sequences have been obtained and tagged from [RFAM](https://rfam.org).

## Alignments 
Each ribosomal RNA has an associated MSA. The alignments were generated according to the procedure described at (https://doi.org/10.1093/molbev/msy101)

## 2D maps 
Pregenerated 2D layouts of RNAs were exported into the RNA topology viewer using the [EMBL-EBI PDBe API](https://www.ebi.ac.uk/pdbe/static/entry/7k00_1_A.json).

## 3D Structures 
3D structures were fetched from the RCSB using the APIs of [RCSB coordinate server](https://models.rcsb.org/v1/7k00/).

## 3D coloring
To facilitate coloring of the 3D structures, the RiboVision 2.0 color themes were created within PDBe Mol*, using the color wrapper of (Proteopedia)[https://github.com/molstar/molstar/tree/master/src/examples/proteopedia-wrapper]  as a template. The tailored Mol* code is available from local [GitHub repo](https://github.com/LDWLab/pdbe-molstar-GT).

## Sequence and structure associated data
### Protein contacts
RNA-Protein contacts are computed upon selecting a specific RNA complex and specifying the main RNA chain. The contacts are computed by [NeighborSearch module of BioPython](https://biopython.org/docs/1.75/api/Bio.PDB.NeighborSearch.html) using the [KD Tree algorithm](https://biopython.org/docs/1.75/api/Bio.KDTree.KDTree.html) with a cutoff distance of 3.5 A.

### Chemical modifications
Chemical modifications (if any) are extracted from "_entity_poly.pdbx_seq_one_letter_code" fields of the selected CIF file. 

## Available attributes for calculated mapping data:

### Shannon Entropy
The Shannon entropy (as well as all properties listed below) was computed from the gap adjusted probabilities as:

<img src="https://render.githubusercontent.com/render/math?math=H_{SE}(n) = -\sum_{i=1}^c p_i(n)log_2p_i(n) \approx -\sum_{i=1}^c f_i(n)log_2f_i(n)">

### Two group comparison (TwinCons)
In the case of two groups being selected in the phylogeny browser, RiboVision2.0 provides an additional option to compute an in-house developed score, TwinCons. TwinCons (
https://doi.org/10.1371/journal.pcbi.1009541) is computed for a single position of the MSA and compares two pre-defined groups (represented by vectors of the gap adjusted nucleotide frequencies) based on their similarity defined by the pre-computed substitution matrix, blastn (https://github.com/LDWLab/TwinCons). TwinCons represents the transformation price between the two vector columns related by the substitution matrix. 

### Color Schemes
Each calculated attribute is mapped on a matplotlib colorscheme. For single continuum attributes (e.g. Shannon entropy), RiboVision 2.0 uses single continuum colormaps (viridis). For diverging data attributes (e.g. TwinCons), RiboVision 2.0 uses diverging colormaps like Blue-White-Red. All colormaps were generated with the python matplotlib library and exported to JavaScript with the [js-colormaps package](https://github.com/timothygebhard/js-colormaps). Further information about colormaps in [matplotlib](https://bids.github.io/colormap/).
