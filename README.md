# Taxconverter
A lightweight tool for one purpose only: to unify the outputs of different taxonomic classifiers. Currently supports MMseqs2, Centrifuge v1.0.4, Kraken2 v2.1.3, Metabuli v1.1.0, MetaMaps v.633d2e. The output files of these tools are converted to `contigs    predictions` format.

The tab-delimitered output file has two columns: sequence identifiers (the 1st column, `contigs`) and taxonomic labels on all levels from domain to species concatenated with ";" (`predictions`). 

IMPORTANT: if you come from the Taxometer README page (https://github.com/RasmussenLab/vamb/blob/taxometer_release/README_Taxometer.md), run the command with the `--mmseqs-format` flag. The output file will have a different format but the same information. Taxometer, when run from the release branch (`taxometer_release`), uses this format. It is also possible to run Taxometer from the newest release of VAMB, in which case you don't need the flag. This compatibility issue will be fixed in the future releases of the VAMB library.

Suggestions and contributions are most welcome.

## Installation

The package is compatible with Python version <=3.11.

```
pip install taxconverter
```

Or clone this repo and install the package from the source. 

```
git clone git@github.com:RasmussenLab/taxconverter.git
cd taxconverter
pip install -e .
```

## Usage
To convert Centrifuge, Kraken2, MetaMaps and MMSeqs2 outputs, provide one file with the taxonomy annotation results:

```
taxconverter centrifuge -i centrifuge_annotations.tsv -o result.tsv
```

```
taxconverter kraken2 -i kraken2_annotations.tsv -o result.tsv
```

```
taxconverter metamaps -i metamaps_annotations.tsv -o result.tsv
```

```
taxconverter mmseqs2 -i metamaps_annotations.tsv -o result.tsv
```

To convert a Metabuli output, provide two files with `_classifications.tsv` and `_report.tsv` postfixes:

```
taxconverter metabuli -c metabuli_classifications.tsv -r metabuli_report.tsv -o result.tsv
```

For more help, run `taxconverter -h`, `taxconverter metabuli -h`, `taxconverter centrifuge -h`, `taxconverter kraken2 -h`, `taxconverter metamaps -h`, `taxconverter mmseqs2 -h`

## References and links
This package is made to complement the Taxometer tool for refining taxonomic annotations from any classifier using contigs k-mers and co-abundances ([Article](https://www.nature.com/articles/s41467-024-52771-y) [Youtube](https://youtu.be/9vuMs-n1-yU?si=3LD2ayyhET1BeqnC)). 

Other links:
* __MMseqs2__ [Article](https://academic.oup.com/bioinformatics/article/37/18/3029/6178277?login=true) [Github](https://github.com/soedinglab/MMseqs2)
* __Metabuli__ [Article](https://www.biorxiv.org/content/10.1101/2023.05.31.543018v2) [Github](https://github.com/steineggerlab/Metabuli)
* __Centrifuge__ [Article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5131823/) [Github](https://github.com/infphilo/centrifuge)
* __Kraken2__ [Article](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) [Github](https://github.com/DerrickWood/kraken2)
* __MetaMaps__ [Article](https://www.nature.com/articles/s41467-019-10934-2) [Github](https://github.com/DiltheyLab/MetaMaps)