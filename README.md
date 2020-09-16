# GCsnap

GCsnap is a flexible Python-based tool that allows for the interactive comparison of the genomic contexts of protein-coding genes from any genome at any taxonomic level, integrating them with functional and structural information for any of the genes shown. 

![figure1](https://github.com/JoanaMPereira/GCsnap/blob/master/examples/Fig1.png)

By connecting the output to different protein databases, the user can navigate through the different genomic contexts from a simple interactive platform, facilitating the further analysis of the contexts found. 

GCsnap is not limited to a single input format, can preform batch jobs and accepts protein classification maps. 

All information is stored in detailed, human and machine-readable files, and customable publication-ready figures.

## Dependencies

GCsnap is written in Python 3.7.4 and requires mostly core Python modules. Only three external packages are required: 
  - Biopython
  - Bokeh
  - Networkx 

Additionally, GCsnap relies on a local installation of BLASTp and PsiBlast (versions 2.4.0+ and above). 

## Allowed inputs

GCsnap takes as main input a list of sequence identifiers, which can be in **Entrez, UniprotKB, and UniRef formats, or a mix**. These identifiers can be given as:
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
  -n_flanking: the number of flanking genes (to each side) to be taken. By default, it is set to 4.
  -n_flanking5: the number of flanking genes to be taken on the 5' side. By default, it is set to 4.
  -n_flanking3: the number of flanking genes to be taken on the 3' side. By default, it is set to 4.
  -get_taxonomy: set to false if no taxonomy is to be collected.
  -annotate_TM: set to true to annotate the presence of transmembrane segments and signal peptides.
  -annotation_TM_mode: the mode to use to collect transmembrane and signal peptide annotations (phobius, tmhmm or uniprot).
  -clans_pattern: a set of patterns in CLANS groups names that define different groups to be considered as independent jobs.
```
### 1. Simple job

Using the example in folder `example/ybez_KHI`, the input file `targets_ybez_selected.txt` contains a list of protein sequence identifiers in UniprotKB format. Running:
```
python3 GCsnap.py -targets targets_ybez_selected.txt
```
will generate the output folder `targets_ybez_selected`, where all output files and figures are stored.
This will NOT annotate transmembrane segments and signal peptides.

In order to do so, one shall run:
```
python3 GCsnap.py -targets targets_ybez_selected.txt -annotate_TM True
```
which will by default collect that information from Uniprot.

### 2. Job from a CLANS file

Using the example in folder `examples/yqlc_KHII/`, the input file `yqlc_nostoc_blast_nrbac70.clans` is a classification file encompassing two clusters of sequences, which are named `cluster1_cyanobacteria` and `cluster2_allothers`. 
Running:
```
python3 GCsnap.py -targets yqlc_nostoc_blast_nrbac70.clans 
```
will make GCsnap consider all identifiers as a single job, while running:
```
python3 GCsnap.py -targets yqlc_nostoc_blast_nrbac70.clans -clans_pattern cluster
```
will make GCsnap identify all clusters in the CLANS file that have 'cluster' in their name, which will be considered as two independent jobs, generating the two folders `cluster1_cyanobacteria` and `cluster2_allothers`.

The figure below depicts the cluster map and the 'most conserved genomic feature' panels from each of these jobs.
![figure2](https://github.com/JoanaMPereira/GCsnap/blob/master/examples/Fig2.png?v=4&s=200)

## Acknowledgements

GCsnap was developed during the COVID-19 lockdown. 

I would like to thank Prof. Andrei N. Lupas, Dr. Laura Weidmann-Krebs, Dr. Marcus Hartmann, Dr. Vikram Alva, Dr. Felipe Merino, Dr. Jörg Martin, Dr. Adrian Fuchs, Hadeer Elhabashy, Dr. João Barros and Prof. Volkmar Braun for the great support and the insightful discussions that helped the development of GCsnap.
