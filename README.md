# RNA_virus_detector. 
This is an RNA virus identification tool. Currently in preparation.  

Dependency. 
  Prodigal  
  Hmmer  
  python3  
    
Usage:  
  The predict.py does not accept any command line input now. For testing, please modify the file path at the end of the code    
    filepath : your data work space(test)  
    filename1 : the contigs file name (Ex: test.fa)  
    filename2 : the hmmer search result (Ex: hmmer_search)  
    filename3 : the predicted proteins file name (Ex: test.faa)  
    filename4 : the predicted positive contigs file (Ex: pos_contigs.vb.fasta)  
    
  For testing:  
```
    cd test. 
    prodigal -p meta -i test.fa -a test.faa. 
    hmmsearch --tblout hmmer_search --noali -E 0.001 --cpu 112 ../ref/RNA_virus.hmm test.faa. 
    python3 predict.py. 
```
