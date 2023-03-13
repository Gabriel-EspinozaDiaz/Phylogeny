# Phylogeny
Uses 'Unweighted Pair Group Method with Arithmetic Mean' (UPGMA) to deduce how closely related species are based on the amino acid sequences of ancient genes.

This project is in it's very basic stages at the moment, as at the moment there is no actual phylogenetic tree being deduced, just an illustration of the similiarity matrix used to make the tree, with each iteration stored as a numpy matrix and printed at each repetition of the UPGMA cycle. 

The current goal now is to rearrange the program into a class-based system. The subsequent plan to finish the program goes as follows:
1. Store all saved values with their respective cluster
2. Store data excised from the table in some sort of tree data structure
3. Illustrate the diagram (this most likely will be done with pyplot)
