# Taxconverter
A lightweight tool for one purpose only: to unify the outputs of different taxonomic classifiers. Currently supports Centrifuge v1.0.4, Kraken2 v2.1.3 and Metabuli v1.0.1. The output files of these tools are converted to MMseqs2 format. 

Having everything converted to MMseqs2 format means the output file will only have two informative columns: sequence identifiers (the 1st column) and taxonomic labels on all levels from domain to species concatenated with ";" (9th column). The rest is filled with zeros or by trivial parsing. 

Suggestions and contributions are most welcome.

## Installation
1) Clone this repo and install the package from the source (releasing `pip` package WIP).

```
git clone git@github.com:RasmussenLab/taxconverter.git
cd taxconverter
pip install -e .
```

2) Unzip the two files from `data/lineage.zip` (38.3 MB): `ncbi_lineage.csv` (246.2 MB) and `metabuli_lineage.csv` (58.1 MB), and place them to the `data/` folder.

## Usage
To convert Centrifuge and Kraken2 outputs, provide one file with the taxonomy annotation results:

```
taxconverter centrifuge -i centrifuge_annotations.tsv -o result.tsv
```

```
taxconverter kraken2 -i kraken2_annotations.tsv -o result.tsv
```

To convert a Metabuli output, provide two files with `_classifications.tsv` and `_report.tsv` postfixes:

```
taxconverter metabuli -c metabuli_classifications.tsv -r metabuli_report.tsv -o result.tsv
```

For more help, run `taxconverter -h`, `taxconverter metabuli -h`, `taxconverter centrifuge -h`, `taxconverter kraken2 -h`

## References and links
This package is made to complement the Taxometer tool for refining taxonomic annotations from any classifier using contigs k-mers and co-abundances (link). 

Other links:
* __MMseqs2__ [Article](https://academic.oup.com/bioinformatics/article/37/18/3029/6178277?login=true) [Github](https://github.com/soedinglab/MMseqs2)
* __Metabuli__ [Article](https://www.biorxiv.org/content/10.1101/2023.05.31.543018v2) [Github](https://github.com/steineggerlab/Metabuli)
* __Centrifuge__ [Article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5131823/) [Github](https://github.com/infphilo/centrifuge)
* __Kraken2__ [Article](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) [Github](https://github.com/DerrickWood/kraken2)