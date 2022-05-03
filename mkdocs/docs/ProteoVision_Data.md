# ProteoVision Data
## Phylogeny (DESIRE)
The subset of 179 species from the [DESIRE](https://doi.org/10.1093/molbev/msy101) (Sparse and Efficient Representation of Extant Biology), 
database was organized into a phylogenetic browser using a tree topology from [NCBI](https://www.ncbi.nlm.nih.gov/taxonomy).

## Alignments 
Each ribosomal protein has an associated MSA. First, an MSA reference was generated with [MATRAS](https://doi.org/10.1093/nar/gkg581) from multiple structure superimpositions. Then, amino acid sequences of species from the DESIRE database were added to the reference alignment using [MAFFT](https://doi.org/10.1093/bioinformatics/bts578).

## 2D maps 
Topologies of the protein secondary structures (Laskowski; [10.1093/nar/gkn860](https://dx.doi.org/10.1093%2Fnar%2Fgkn860)) were exported into PDB topology viewer using the [EMBL-EBI PDBe API](https://www.ebi.ac.uk/pdbe/api/doc/).

## 3D Structures 
3D structures were fetched from the PDBe using the APIs of [EMBL-EBI coordinate server](https://www.ebi.ac.uk/pdbe/coordinates/). The selection of ranges was implemented using the syntax of the [LiteMol’s coordinate server](https://coords.litemol.org/).

## Alignment associated data (Fold, Phase)
DESIRE holds annotations of domain architecture from ECOD (Cheng; [10.1371/journal.pcbi.1003926](https://doi.org/10.1371/journal.pcbi.1003926)) and ribosomal phase definitions (Kovacs; [10.1093/molbev/msx086](https://doi.org/10.1093/molbev/msx086)) at residue level for one representative species (*E. coli*). Using the alignments, ProteoVision retrieves these annotations for each column of an alignment and displays them as a hovering pop-up next to each residue.

## Available attributes for calculated mapping data:
### Amino Acid frequencies
Amino acid frequencies in each column of an MSA were adjusted for presence of gaps. Thus, the gap frequencies were prorated and were treated as a uniform distribution among all possible amino acid characters, such that a single character in a gap counts as 0.05, as described by [Bernier et al.](https://doi.org/10.1093/molbev/msy101).

### Shannon Entropy
The Shannon entropy (as well as all properties listed below) was computed from the gap adjusted probabilities as:

<img src="https://render.githubusercontent.com/render/math?math=H_{SE}(n) = -\sum_{i=1}^c p_i(n)log_2p_i(n) \approx -\sum_{i=1}^c f_i(n)log_2f_i(n)">

### Two group comparison (TwinCons)
In case of two groups selected in the phylogeny browser, ProteoVision provides an additional option to compute an in house developed score, TwinCons. TwinCons is computed for a single position of the MSA that compares two pre-defined groups (represented by vectors of the gap adjusted amino acid frequencies) based on their similarity defined by the pre-computed substitution matrix. TwinCons represents the transformation price between the two vector columns related by the substitution matrix.

### Charge, hydropathy, hydrophobicity, polarity, mutability
The physico-chemical properties for each position within an MSA are computed as average properties for a given distribution of the amino acid frequencies. The tabulated values for each property were obtained from the available literature:
- [charges](https://doi.org/10.1186/1758-2946-5-39)
- [hydropathy](https://doi.org/10.1016/0022-2836(82)90515-0)
- [hydrophobicity](https://doi.org/10.1093/protein/5.5.373)
- [polarity](https://doi.org/10.1016/0022-5193(68)90069-6)
- [mutability](https://doi.org/10.1093/bioinformatics/8.3.275)

### Color Schemes
Each calculated attribute is mapped on a matplotlib colorscheme. For single continuum attributes (like Shannon entropy or Polarity), ProteoVision uses single continuum colormaps like plasma and viridis. For diverging data attributes (like Charge or TwinCons), ProteoVision uses diverging colormaps like Blue-White-Red or Green-White-Purple. All colormaps were generated with the python matplotlib library and exported to JavaScript with the [js-colormaps package](https://github.com/timothygebhard/js-colormaps). Further information about colormaps in [matplotlib](https://bids.github.io/colormap/).
