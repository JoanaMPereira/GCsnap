# GCsnap

GCsnap is a flexible Python-based tool that allows for the interactive comparison of the genomic contexts of protein-coding genes from any genome at any taxonomic level, integrating them with functional and structural information for any of the genes shown. 

![figure1](https://github.com/JoanaMPereira/GCsnap/blob/master/examples/Fig1.png)
> Figure 1. Interactive output for the comparison of the genomic contexts of bacterial ribosomal protein S16 (rpsP). (A) The core panel illustrating the different genomic contexts centred in the target inputs (dark grey). Genes are represented as arrow boxes and coloured based on the family clusters identified. When hovering over a gene box, an information box that summarises the information collected for that given gene appears. (B) Legend for the families shown in A. When hovering over a gene box or its description, an information box summarising the structural and functional information collected for each member of that family appears. Clicking in a box redirects the user to the page of a representative from that family in the SWISS-MODEL repository. (C) Cladogram representing the hierarchy used to sort the labels in A. The tips store the taxonomic information collected for the genomic context represented. (D) Summarising tabs. Three different tabs are shown, which summarise the genomic contexts shown in three different ways. (E) Figure tools from Bokeh. These buttons allow the user to activate (or deactivate) specific interactive features (e.g., activate zooming into specific regions) and to save each panel as individual figures. 

By connecting the output to different protein databases, the user can navigate through the different genomic contexts from a simple interactive platform, facilitating the further analysis of the contexts found. GCsnap is not limited to a single input format, can preform batch jobs and accepts protein classification maps. 

All information is stored in detailed, human and machine-readable files, and customable publication-ready figures.

## Dependencies

GCsnap is written in Python 3.7.4 and requires mostly core Python modules. Only three external packages are required: 
  - Biopython
  - Bokeh
  - Networkx 

Additionally, GCsnap relies on a local installation of BLASTp and PsiBlast (versions 2.4.0+ and above). 
