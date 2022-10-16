# GCsnap

GCsnap is a flexible Python-based tool that allows for the interactive comparison of the genomic contexts of protein-coding genes from any genome at any taxonomic level, integrating them with functional and structural information for any of the genes shown. 

![figure1](https://github.com/JoanaMPereira/GCsnap/blob/master/examples/Fig1.png)

By connecting the output to different protein databases, the user can navigate through the different genomic contexts from a simple interactive platform, facilitating the further analysis of the contexts found. 

GCsnap is not limited to a single input format, can preform batch jobs and accepts protein classification maps. 

All information is stored in detailed, human and machine-readable files, and customable publication-ready figures.

## Dependencies

GCsnap is written in Python 3.7 and should run on Python 3.x. It was tested on Python 3.7 and 3.8. It requires mostly core Python modules and only five external packages are required: 
  - Biopython
  - Bokeh
  - Networkx 
  - PaCMAP
  - Scikit-learn

For detailed requirements, check ```requirements.txt```.

Additionally, GCsnap relies on a local installation of BLASTp and PsiBlast (versions 2.4.0+ and above) and MMseqs. 

## Installation

### Installing from Source

Download the zip archive or clone the repository with git:

```
# To download
git clone https://github.com/JoanaMPereira/GCsnap
cd GCsnap

# To update
git pull origin master

# To install
python setup.py install
```

## Allowed inputs

GCsnap takes as main input a list of sequence identifiers, which can be in **Entrez, UniprotKB, UniRef, GeneID, and ENSEMBLE ID formats, or a mix**. These identifiers can be given as:
  - a text file, where each is in a different line
  - a fasta file, where the sequence header starts with the sequence identifier
  - a sequences cluster file in CLANS format
  - direct input in the terminal as a space-separated list
  
## Usage

In its most simple mode of usage, GCsnap only requires a list of sequence identifiers. 

**Required** arguments are:
```
  -targets: which can be a list of sequence identifiers, a text file, a fasta file, a clans file or a list of files
```
**Optional** arguments allow for the tweaking of GCsnap behaviour. There are various optional arguments that can be used, but the most relevant are:
```  
  -user_email: it may be required to access the NCBI databases. It is not used for anything else.
  -ncbi_api_key: the key for NCBI API, which allows for up to 10 queries per second to NCBI databases. Can be obtained after obtaing an NCBI account.
  -cpu: the number of cpus used for running. By default, it is set to 1. Using more allows GCsnap to parallelize the most time-consuming steps. We recommend using this in combination with a ncbi api key.
  -n_flanking: the number of flanking genes (to each side) to be taken. By default, it is set to 4.
  -n_flanking5: the number of flanking genes to be taken on the 5' side. By default, it is set to 4.
  -n_flanking3: the number of flanking genes to be taken on the 3' side. By default, it is set to 4.
  -get_taxonomy: set to false if no taxonomy is to be collected.
  -annotate_TM: set to true to annotate the presence of transmembrane segments and signal peptides.
  -annotation_TM_mode: the mode to use to collect transmembrane and signal peptide annotations (phobius, tmhmm or uniprot).
  -clans_pattern: a set of patterns in CLANS groups names that define different groups to be considered as independent jobs.
  -operon_cluster_advanced: set to true to have a more comprehensive analysis/summary of the genomic contexts found. Ideal for very large input sets.
  -all-against-all_method: BLASTp is as default, but can be set to mmseqs if MMseqs2 is available.
```
### 1. Simple job

Using the example in folder `example/ybez_KHI`, the input file `targets_ybez_selected.txt` contains a list of protein sequence identifiers in UniprotKB format. Running:
```
GCsnap -targets targets_ybez_selected.txt
```
will generate the output folder `targets_ybez_selected`, where all output files and figures are stored.
This will NOT annotate transmembrane segments and signal peptides.

In order to do so, one shall run:
```
GCsnap -targets targets_ybez_selected.txt -annotate_TM True
```
which will by default collect that information from Uniprot.

### 2. Job from a CLANS file

Using the example in folder `examples/yqlc_KHII/`, the input file `yqlc_nostoc_blast_nrbac70.clans` is a classification file encompassing two clusters of sequences, which are named `cluster1_cyanobacteria` and `cluster2_allothers`. 
Running:
```
GCsnap -targets yqlc_nostoc_blast_nrbac70.clans 
```
will make GCsnap consider all identifiers as a single job, while running:
```
GCsnap -targets yqlc_nostoc_blast_nrbac70.clans -clans_pattern cluster
```
will make GCsnap identify all clusters in the CLANS file that have 'cluster' in their name, which will be considered as two independent jobs, generating the two folders `cluster1_cyanobacteria` and `cluster2_allothers`.

The figure below depicts the cluster map and the 'most conserved genomic feature' panels from each of these jobs.
![figure2](https://github.com/JoanaMPereira/GCsnap/blob/master/examples/Fig2.png)

### NEW: 3. Advanced genomic context analysis

Since version 1.0.16, GCsnap incorporates an "Advanced mode", which can be used by setting the `-operon_cluster_advanced` flag to True. This mode uses `PacMap` (https://github.com/YingfanWang/PaCMAP) to identify clusters of similar genomic contexts and, instead of displaying all contexts in one single view, first generates a summary page displaying the clusters of genomic contexts found as well as a "family composition spectrum" that allows for an easier interpretation of the diversity of genomic contexts identified. In addition, a detailed page is generated for each genomic context type defined, which includes the classic genomic context block but also a table listing the properties of families found. This is useful for very large input sets (thousands of sequences or large clans maps)

Using the example in folder `examples/yqlc_KHII/` and the input file `yqlc_nostoc_blast_nrbac70.clans` without defining target clusters (i.e., it will use all sequences in the map), we can run a simple advanced job by running:

```
GCsnap -targets yqlc_nostoc_blast_nrbac70.clans -get_taxonomy True -operon_cluster_advanced True -get_pdb True -get_functional_annotations True -interactive True
```

The output summary page displays on the top 3 different scatter plots, where the color of the dots corresponds to the type of the genomic context they belong to. As the output file was a CLANS file, GCsnap recognized it as such and so diplays it. Otherwise, the area would be white unless an input CLANS file is given explicitly.

![figure3](https://github.com/JoanaMPereira/GCsnap/blob/master/examples/Fig3.png)

Each sequence is listed below in a table, with their corresponding associated genomic context type as well as the collected taxonomy. The button next to the top of the table allows for changing the display to show a sorted spectrum of family compositions for each context type:

The figure below shows the summary page for genomic context type 0000.
![figure4](https://github.com/JoanaMPereira/GCsnap/blob/master/examples/Fig4.png)

## Citing GCsnap

GCsnap is going to be published in the Computation Resources (2021) Special Issue of Journal of Molecular Biology (JMB) and you can cite it already:

```
J. Pereira, GCsnap: interactive snapshots for the comparison of protein-coding genomic contexts, 
J. Mol. Biol. (2021) 166943. https://doi.org/https://doi.org/10.1016/j.jmb.2021.166943.
```

## Acknowledgements

GCsnap was developed during the COVID-19 lockdown. 

I would like to thank Prof. Andrei N. Lupas, Dr. Laura Weidmann-Krebs, Dr. Marcus Hartmann, Dr. Vikram Alva, Dr. Felipe Merino, Dr. Jörg Martin, Dr. Adrian Fuchs, Hadeer Elhabashy, Prof. Volkmar Braun, Dr. João Rodrigues and Dr. João Barros for the great support and the insightful discussions that helped the development of GCsnap.

GCsnap is being maintained at the Biozentrum of the University of Basel, and I would like to thank the Schwede team and the Basler group for insightful discussions that are driving many new developments.
