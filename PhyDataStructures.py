import numpy as np
import math

'''
- Reads in data 
- Constructs a table using unweighted pairwise distance by methods of arithmetic means (UPGMA)
- 
'''
class Pairwise:
    def __init__(self,directory):
        if type(directory) != str:
            raise TypeError('DIRECTORY MUST BE IN STRING FORM')
        phy_file = open(directory,'r')
        phy_string = phy_file.read()
        phyList = phy_string.split('\n')
        
        self.specDict = {}
        self.specList = []

        for n in range(0,len(phyList),2): #Generates a dictionary with species names as keys and sequences as items
            self.specList.append(phyList[n])
            self.specDict.update({phyList[n]:phyList[n+1]})

        specNum = len(self.specDict)
        
        self.parent_table = np.zeros((specNum,specNum),dtype=float)
        self.grouped_table = self.parent_table.copy()

        self.groups = np.array([])

    def checkSequences(sequences):
        '''
        Ensures that all inputs are correct
        Raises any errors present
        '''
        allowed_AAs = ['A','R','N','D','C','Q','E','G','H',
    'I','L','K','M','F','P','O','S','U','T','W','Y','V',]
        for n in range(len(sequences)): #Ensures that all AA sequences are the same length
            if len(sequences[n]) != len(sequences[0]):
                raise TypeError('ALL SEQUENCES MUST BE SAME LENGTH')
    
        for n in range(len(sequences)): #Ensures that AA sequences don't contain any incorrect AA 1-letter codes
            for letter in range(len(sequences[0])):
                if sequences[n][letter] not in allowed_AAs:
                    raise TypeError('SEQUENCES(S) CONTAIN INVALID AMINO ACID CODES')

    def differences(self,spec1,spec2):
        '''
        Takes the index of 2 species and calculates and returns the number of pairwise differences.
        Takes into account 
        '''
        spec1 = self.specDict[self.specList[spec1]]
        spec2 = self.specDict[self.specList[spec2]]
        count = 0
        for n in range(len(spec1)):
            if spec1[n] != spec2[n]:
                count += 1
        return round(math.log(1-(count/len(spec1)))*-100,2)
    
    def constructTable(self):
        '''
        builds the original table
        '''
        table = self.parent_table
        exc = 0
        for n in range(len(self.specList)):
            exc += 1
            for i in range(exc,len(self.specList)):
                table[n][i] = self.differences(n,i)
        return table

    def findMin(self):
        '''
        Runs through the table and finds the smallest difference
        returns the array's indices for minimum difference
        '''

        return 0


'''

'''
class PhyTree:

    def __init__(self):
        pass
