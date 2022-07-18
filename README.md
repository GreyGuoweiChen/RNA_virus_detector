<<<<<<< HEAD
<<<<<<< HEAD
# RNA_virus_detector  
This is an RNA virus identification tool. Currently in preparation.  

## Dependency:
  Prodigal  
  Hmmer  
  python3  
    
## Usage:
  The predict.py does not accept any command line input now. For testing, please modify the file path at the end of the code    
 >   filepath : your data work space(Ex: test)  
 >   filename1 : the contigs file name (Ex: test.fa)  
 >   filename2 : the hmmer search result (Ex: hmmer_search)  
 >   filename3 : the predicted proteins file name (Ex: test.faa)  
 >   filename4 : the predicted positive contigs file (Ex: pos_contigs.vb.fasta)  
    
  For testing:  
```
    cd test  
    prodigal -p meta -i test.fa -a test.faa  
    hmmsearch --tblout hmmer_search --noali -E 0.001 --cpu 112 ../ref/RNA_virus.hmm test.faa  
    python3 predict.py  
=======
=======
>>>>>>> 1a0edebce622e9594af96c685b6c4b1e349b6777
# VirBot: a protein-based RNA virus detector for metagenomic data
VirBot is an RNA virus detection tool, which allows accurate and sensitive identificaton for RNA viral contigs. Currently still in development, feel free to use and welcome any suggestions.  

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
<<<<<<< HEAD
>>>>>>> 8482fc5 (submit for testing 2)
=======
>>>>>>> 1a0edebce622e9594af96c685b6c4b1e349b6777
```
