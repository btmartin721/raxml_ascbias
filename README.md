# extractVariantSites4raxml  
Python 3 script for removing invariant sites from a PHYLIP file to use the RAxML 8.2.X ascertainment bias correction  

The script removes sites that are invariant. RAxML considers IUPAC characters and gaps to be invariant if phasing them could yield invariant sites. Consider the following example:  

```  
AAT  
ART  
TAG  
T-G  
```  

extractVariantSites.py will remove column 2 and yield:  

```  
AT  
AT  
TG  
TG  
```  

Usage: `./extractVariantSites.py -f <INFILE> -o <OUTFILE; default="out.phy">`  

The program has three python3 dependencies: numpy, pandas, and biopython  
You can install them using Anaconda or as follows:  

`sudo apt-get update`  
`sudo apt-get install python3-numpy`   
`sudo apt-get install python3-pandas`  
`sudo apt-get install python3-biopython`  

They can also be installed using pip:  

`sudo apt-get install python3-pip`  
`sudo pip3 install pandas`  

etc.  




