FluorosequencingSimulation
==========================



INTRODUCTION
____________________________

The scripts (a) seqmodlibMP.py and (b) randsiggen.c codes for the steps in the Monte Carlo simulation detailed in the paper. The two scripts form the underlying code and a wrapper code must be applied according to the parallelized computing infrastructure used. 

System Requirements
____________________________

1. Python 2.7 
2. Standard C compiler 
3. High Performance Computing (HPC) infrastructure capable of parallel computing.
4. libGSL
5. SMFT Merseinne Twister

INSTALLATION
____________________________

Unzipr or Untar the files into a folder (preferably added to the PYTHONPATH)
*Ipython IDE may be installed to write and execute the wrapper codes
randsiggen.c needs to know where libGSL and SMFT header files are


PROGRAMME EXECUTION
____________________________

(a) PREPROCESSING PROTEOME FILE

	1. Proteome of choice can be downloaded from a number of publicly available databases like Uniprot (http://www.uniprot.org/) etc user defined. 
	2. The file is parsed to yield a python dictionary in the form : 
		 {
	           'PROTEIN NAME 1': 'AMINO ACID SEQUENCE 1',
	           'PROTEIN NAME 2': 'AMINO ACID SEQUENCE 2',
        	 ...}
	3. The file is stored as a pkl file. 

(b) WRAPPER CODE

	For example: 
	
		import seqmodlibMP
		proteome = seqmodlibMP.load_proteome(proteome_pkl_File) #proteome pkl file is loaded
		cleaved = seqmodlibMP.cleave(proteome,'E') # This cleaves all peptides after 'E' (representing the GluC cleavage)
		attached = seqmodlibMP.attach(cleaved,'C') # This retains only 'C' containing peptides 
		
		trie = seqmodlibMP.monte_carlo_trie(attached, 
						p=0.9,    #90% Edman efficiency
						b=0.0033, #photobleaching rate constant of 0.0033
						u=0.20,   #fluorophore failure rate of 20%
						windows={'K':range(1,30)}, #Considers upto 30 experimental cycles with the label on 'K'
						sample_size = 1000, #Represents the simulation depth
						random_seed=random.random(), #Generates the appropriate random model for the error sources
						silent=True #True turns off the print statement. Useful for debugging
						    )

		trie.find_uniques() 	# Trie Analysis is done by calling 

**NOTE**
*This module also needs to be implemented in a HPC (see system requirements). The function needs a wrapper that partitions the resultant trie over a number of multiterabyte-sized tries. *

(c) RESULT FILE

	The end result file has a mapping of a fluorosequence to the tally of source peptides generating it
	For example: 
		xxKxxk : (Protein A, 20),(Protein B, 8), ...
		xKxxK  : (Protein C, 200), (Protein A, 20), ...
		...

(d) VISUALIZING DATA

All the data were graphed for visualization using standard python matplotlib library functions (http://matplotlib.org/)



AUTHOR
_____________________________
Alexander Boulgakov (alexander[DOT]boulgakov{AT}utexas[DOT]edu)




  
