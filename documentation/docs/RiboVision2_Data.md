# RiboVision 2 Data
## Phylogeny (DESIRE)
The subset of 152 species from the [DESIRE](https://doi.org/10.1093/molbev/msy101) (Sparse and Efficient Representation of Extant Biology), 
database was organized into a phylogenetic browser using a tree topology from the [Banfiled lab](https://doi.org/10.1038/nmicrobiol.2016.48).

## Alignments 
Each ribosomal RNA has an associated MSA. The alignments were generated accodring to procedure desribed at (https://doi.org/10.1093/molbev/msy101)

## 2D maps 
Topologies of the protein secondary structures (Laskowski; [10.1093/nar/gkn860](https://dx.doi.org/10.1093%2Fnar%2Fgkn860)) were exported into PDB topology viewer using the [EMBL-EBI PDBe API](https://www.ebi.ac.uk/pdbe/api/doc/).

## 3D Structures 
3D structures were fetched from the PDBe using the APIs of [EMBL-EBI coordinate server](https://www.ebi.ac.uk/pdbe/coordinates/). The selection of ranges was implemented using the syntax of the [LiteMol’s coordinate server](https://coords.litemol.org/).

## Secuence and structure associated data (Chemical modifications, Protein Contacts)
Protein contacts are calculated on the fly from the specified RNA-protein complex using the NeighborSearch algorithm implemented in Biopython. All protein chains that are in contact within a cut off distance of 3.5 A of a specified RNA chain are reported in the Main Naviagtion panel upon data processing. The mapping of RNA-protein contacts onto 2D and 3D rperesentation of RNA molecules are performed by selecting the (sub-) set of preprocessed Protein chains from the Main Navigation papnel. Upon selection, RNA residues that are in contact with a given protein chain will be highlighed in a separate color in 2D and 3D templates. Arritionally, a more detailed information regarding the protein contact (protein name and the chain ID) is available via interactive Tooplip implemented in the 2D RNA viewer; the selected proteins will also appear as additional 3D objects in the 3D viewer and colored according to their contact map colors.

Chemical modification  are parsed from  _entity_poly.pdbx_seq_one_letter_code of the CIF file. 
.

## Available attributes for calculated mapping data:
### Nucleotide frequencies
Nucleotide frequencies in each column of an MSA were adjusted for presence of gaps. Thus, the gap frequencies were prorated and were treated as a uniform distribution among all possible nucleotide characters, such that a single character in a gap counts as 0.05, as described by [Bernier et al.](10.1093/molbev/msy101).

### Shannon Entropy
The Shannon entropy (as well as all properties listed below) was computed from the gap adjusted probabilities as:

<img src="https://render.githubusercontent.com/render/math?math=H_{SE}(n) = -\sum_{i=1}^c p_i(n)log_2p_i(n) \approx -\sum_{i=1}^c f_i(n)log_2f_i(n)">

### Two group comparison (TwinCons)
In case of two groups selected in the phylogeny browser, RiiboVision2.0 provides an additional option to compute an in house developed score, TwinCons. TwinCons (
https://doi.org/10.1371/journal.pcbi.1009541) is computed for a single position of the MSA that compares two pre-defined groups (represented by vectors of the gap adjusted nucleotide frequencies) based on their similarity defined by the pre-computed substitution matrix, blastn (https://github.com/LDWLab/TwinCons). TwinCons represents the transformation price between the two vector columns related by the substitution matrix. 

### Color Schemes
Each calculated attribute is mapped on a matplotlib colorscheme. For single continuum attributes (like Shannon entropy), RiboVision2.0 uses single continuum colormaps like plasma and viridis. For diverging data attributes (likeTwinCons), RiboVision2.0 uses diverging colormaps like Blue-White-Red. All colormaps were generated with the python matplotlib library and exported to JavaScript with the [js-colormaps package](https://github.com/timothygebhard/js-colormaps). Further information about colormaps in [matplotlib](https://bids.github.io/colormap/).
