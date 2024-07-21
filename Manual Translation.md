The text files available at https://reaction-networks.net/networks/ list reaction networks as sequences of numbers. Each row of the text file encodes one network. The first two numbers give the number of reactions and the number of species, in that order. For example, a 3-species 2-reaction network will begin with “2 3”. The remaining numbers describe the network itself. 

If a network has r reactions and s species, the numbers 0 through r-1 are assigned to the reactions, and the numbers r to s-1 are assigned to the species. For example, in a 2-species 2-reaction network, 0 and 1 represent the reactions while 2 and 3 represent the species. 

The remaining numbers in the sequence should be read in pairs. The order of the numbers indicates whether a species is in the reactant or product complex of a given reaction, and the amount of times a given number pair is repeated indicates the molecularity of the relevant species in the relevant complex. If the species index appears before the reaction index in the number pair, that species is in the reactant complex. If the species index appears after the reaction index, that species is in the product complex. Let’s walk through an example:

3 1 3 0 0 3 0 3 1 3 2 3 2 3

The number pair, (3 1), tell us this is a 1-species, 3-reaction network. So, the numbers 0, 1, 2 represent reactions, and 3 represents the single species present in this network. 

Then, we have the following number pairs pertaining to the first reaction (indexed as 0): (3 0) (0 3) (0 3). The species index appears before the reaction index once, and after it twice. So, this reaction is A → 2A. 

Next, the pair (1 3) for the second reaction tells us that the product complex contains one copy of A, and the reactant complex is empty. So, this reaction is 0 → A. 

Finally, the pairs (2 3) (2 3) tell us the reactant complex is empty and the product complex contains two copies of A. The final reaction is 0 → 2A. The entire network is therefore A → 2A, 0 → A, 0 → 2A.
