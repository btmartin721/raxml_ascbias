# extractVariantSites4raxml  

Python 3 script for removing and counting invariant sites from a PHYLIP file to use for the RAxML 8.2.X ascertainment bias correction  

The script removes all invariant sites and generates output files for the Felsenstein and Stamatakis ascertainment bias corrections. It now does all three steps in a single execution.  

The script is optimized and way, way faster now. I ran it on a PHYLIP file with 107 individuals and 1,094,776 sites in 273 seconds.  **NOTE: It uses a good bit of memory. On my PC it used ~3 or 4 GB of RAM when using a file with >1,000,000 sites.  
 
Felsenstein Correction: Counts number of invariant sites and writes the count to file for input into RAxML  
Stamatkis Correction: Counts invariant sites for each base (A C G T) and writes counts to space-delimited file.  

RAxML considers IUPAC characters and gaps to be invariant if phasing them could yield invariant sites. Consider the following example:  

```  
AAT  
ART  
TAG  
T-G  
```  

ascbias.py will remove column 2 and yield:  

```  
AT  
AT  
TG  
TG  
```  

Usage: `./ascbias.py -p <PHYLIP_INFILE> -o <OUTFILE; default="out.phy">`  

You might need to call the script with python3 to run it. E.g.:  

`python3 ./ascbias.py -p <PHYLIP_INFILE>`  

##Optional argument:  

```
-o [Specify output file]  
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




