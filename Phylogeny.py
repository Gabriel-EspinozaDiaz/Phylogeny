from PhyDataStructures import Pairwise

class main():
    test1 = Pairwise('Sequences.txt')
    print(test1.parent_table)
    print(test1.merge([1,2]))