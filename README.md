# VirBot: a protein-based RNA virus detector for metagenomic data
VirBot is an RNA virus detection tool, which allows accurate and sensitive identificaton for RNA viral contigs. Currently we are still optimizating the codes. But this version is ready to test/use and we welcome any suggestions.

## Dependency:
* Prodigal
* Hmmer
* python 3.x

### Quick install

We highly recommend using `conda` to install all the dependencies.

```bash
# create the environment and install the dependencies
conda create -n virbot -c bioconda prodigal hmmer python=3
# activate the environment
conda activate virbot
```
Please download VirBot by the "git lfs clone" instread of "git clone".
```
git lfs clone https://github.com/GreyGuoweiChen/RNA_virus_detector
```

## Usage:

```
python predict.py [--input INPUT_CONTIG] [--output PREDICTION] [--temp_dir TEMPORAY FOLDER]
```

### Options 

```
--input: The input contig file in fasta format.
--output: The output fasta file containing all the predicted RNA virus contigs.
--temp_dir: The temporary directory used to hold temporary files (default: temp).
```

  ### Example:  
```
python predict.py --input test/test.fa --output test_results.fna
```
