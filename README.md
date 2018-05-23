# extractVariantSites4raxml  

Python 3 script for removing or counting invariant sites from a PHYLIP file to use the RAxML 8.2.X ascertainment bias correction  
The script does the generates output files for the Lewis, Felsenstein, and Stamatkis ascertainment bias corrections  

Lewis Correction: removes sites that are invariant.  **NOTE:** Lewis takes a while to run for large datasets (e.g., ~70 minutes for 250 individuals X 80,000 sites) 
 
Felsenstein Correction: Counts number of invariant sites and writes count to file for input into RAxML  
Stamatkis Correction: Counts invariant sites consisting of A's, C's, G's, and T's and writes counts to space-delimited file  

RAxML considers IUPAC characters and gaps to be invariant if phasing them could yield invariant sites. Consider the following example:  

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

Usage: `./extractVariantSites.py -f <INFILE> -o <OUTFILE; default="out">`  

You might need to call the script with python3 to run it. E.g.:  

`python3 ./extractVariantSites.py -f <INFILE>`  

##Optional arguments:  

```
-o [Specify outfile prefix]  
-l [Lewis correction; removes all invariant sites]  
-f [Felsenstein correction; counts invariant sites and writes count to file]  
-s [Stamatkis correction; Counts invariant sites for each base. e.g., A C G T]  
```

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




