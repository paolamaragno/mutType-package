# mutType package

Starting from a set of single nucleotide variants in VCF format, the corresponding reference genome and a parameter 'context_length' specified by the 
user, the package function 'mutation_type' determines for each mutation the corresponding mutation type 'UP[REF>ALT]DOWN' such that all mutation types 
have C or T as mutated reference base. The overall length of the mutation type is determined by the context_length parameter.
The package function 'count_table' summarizes the mutation types for the set of mutations into a count table that reports the number of mutations per 
mutation type.
Eventually the package function 'graphical_summary' generates a pdf file showing the barplot visualization of all the mutation types with a frequency 
higher than a threshold specified by the user.

Please visualize the vignette of the package for more information:
browseVignettes('mutType')
