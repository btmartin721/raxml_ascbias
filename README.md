# raxml_ascbias  

Python 3 script for removing and counting invariant sites from a PHYLIP file to use for the RAxML 8.2.X ascertainment bias correction  

The script removes all invariant sites and generates output files for the Felsenstein and Stamatakis ascertainment bias corrections. It now does all three steps in a single execution.  

The script is also optimized and way, way faster now. I ran it on a PHYLIP file with 107 individuals and 1,094,776 sites in 273 seconds.  

**NOTE: It uses a good bit of memory. On my PC it used ~3 or 4 GB of RAM when using a file with >1,000,000 sites.**

## Dependencies:

```
numpy
pandas
biopython
```

You will need the Python 3 versions of each dependency. They can be easily installed using Anaconda or using apt-get or pip.

## Output files:
out.phy: PHYLIP file with invariant sites removed.  
out.phy.felsenstein: Counts number of invariant sites and writes the count to file for input into RAxML (Felsenstein correction):  

`#_invariant_sites`  

out.phy.stamatakis: Counts invariant sites for each base and writes counts to space-delimited file:

`A's C's G's T's`

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

## Usage:
`./ascbias.py -p <PHYLIP_INFILE> -o <OUTFILE; optional (default="out.phy")>`  

You might need to call the script with python3 to run it, depending on your system setup. E.g.:  

`python3 ./ascbias.py -p <PHYLIP_INFILE>`  
