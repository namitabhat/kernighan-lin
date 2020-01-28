# kernighan-lin

Implementation of kernighan lin algorithm for graph partitioning in Matlab.

Input:csv file containing edges from node in the 1st column to node in the 2nd column

Output: .gv file that can be used to visualise partitioning

kl algorithm:
The optimal partition will minimise the cut size (number of edges between the two partitions) by iteratively swapping node pairs between the groups. 

It starts with an even partitioning between nodes into group 1 and group 2.The main inner loop iterates through all nodes in group 1, and picks the node pair that reduces the cut size the most when swapped. A node pair here consists of one node from each group. The selected node pair is swapped, and added to a 'forbidden set' so that these nodes are not touched in later iterations. The cut size after each swap and the states of group 1 and group 2 at each iteration is noted. 

After the main inner loop is completed, the smallest cut size during the iterations is chosen and the corresponding state of group 1 and group 2 is the starting point for the next outer loop. The outer loop is run till there is no improvement in cut size. 

Based on reference from https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/kernlin.pdf and http://cs.baylor.edu/~maurer/partitioning.pdf

