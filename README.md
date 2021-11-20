# Bioinformatics_Protein_Structure

1. Coverage.pl 
- find overlaps in a BED file.
- ```Usage: ./coverage.pl <input bed file name> <output file name> ```

- Input: It takes in input the name of a protein with extension .pdb	
- 1.It reads the corresponding PDB file, extracts from it the coordinates of the C_alpha atoms and creates three arrays x, y and z with such coordinates.			   
- 2.It takes in the coordinates, and calls a subroutine that computes the centroid of the atoms	
- 3.It takes in the coordinates, centroid position, and calls a subroutine that computes the distance of every C_alpha atom of the protein from its centroid							    
- 4.It takes in the distances, and calls a subroutine that prints the histogram of the distances of the atoms from the centroid.
- ```Usage: perl Centroid.pl protein_name.pdb```

2. PPI_network.pl 
- Input: the Protein-Protein Interaction (PPI) graph in sif format. Each line of the file represents an edge.        
- Output:                                                                                                            
- 1.a table (histogram) of the degree frequencies, where each line consists of a value of k and the corresponding 𝑝" 
- 2.the average clustering coefficient AVG_C of the network                                                          
- 3.for all k, the average clustering coefficient AVG_C(k) of all nodes of degree k    
- ```Usage: perl PPI_network.pl ppi.sif```

3. Protein-Protein_Interfaces.pl
- 1.It takes as arguments the C_alpha atoms of the amino acids of each of the two chains and computes the distances Dist_Calpha of all pairs of C_alpha atoms, with one C_alpha in chain A and the other in chain B. It outputs the interfaced amino acids                              
- 2.for each chain computes the fraction of the interface amino acids lying on the secondary structures alpha helices and beta sheets 
- 3.for each interface amino acid of a chain determines the closest interface atom of the same chain to the right in the primary sequence. Then it determines the difference in position of the two amino acids in the sequence.  
- ```Usage: perl Protein-Protein_Interfaces.pl protein_name.pdb, chain1, chain2, threshold(distance threshold between two amino acids)```

4. closeness.pl
- Comput The top 10 nodes (sorted) with the highest closeness centrality value
- ```Usage: perl closeness.pl ppi.sif```
