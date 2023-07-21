from PhyDataStructures import Pairwise

class main():
    test1 = Pairwise('Sequences.txt')
    print(test1.differences(0,1))
    test1.constructTable()
    print(test1.grouped_table)
    print(test1.findMin())
    print(test1.scrapeNames('name & name & name & name'))