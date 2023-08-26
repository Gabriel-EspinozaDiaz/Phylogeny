import numpy as np
import math
import re

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
    
        self.specDict = {} #List of species with associated amino acid sequences
        self.specList = [] #List of species
        for n in range(0,len(phyList),2): #Generates a dictionary with species names as keys and sequences as items
            self.specList.append([phyList[n]])
            self.specDict.update({phyList[n]:phyList[n+1]})

        specNum = len(self.specDict) #number of species
        
        self.parent_table = np.zeros((specNum,specNum),dtype=float) #unchanged table with the original values for reference during merging

        self.groups = np.array([]) #Stores species as groups for sorting out clades

        phy_file.close()

    def check_sequences(sequences):
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
        Takes multiple hits into account using the Poisson Formula
        '''
        spec1 = self.specDict[self.specList[spec1]]
        spec2 = self.specDict[self.specList[spec2]]
        count = 0
        for n in range(len(spec1)):
            if spec1[n] != spec2[n]:
                count += 1
        return round(math.log(1-(count/len(spec1)))*-100,2)
    
    def construct_table(self):
        '''
        builds the original table

        '''
        table = self.parent_table
        exc = 0
        for n in range(len(self.specList)):
            exc += 1
            for i in range(exc,len(self.specList)):
                table[n][i] = self.differences(n,i)
        self.grouped_table = self.parent_table.copy()
        
    def merge(self,coords):
        '''
        Uses the parent table and index to average 2 species/group differences
        Accepts 2 indices for groups on the current table
        Returns the merged value
        '''
        gs1 = coords[0]
        gs2 = coords[1]
        dividend = np.empty(len(gs1)*len(gs2))
        # This code looks important but idk why I wrote it lmao: for n in [i for i in range(len(self.specDict)) if n != gs1 or n != gs2]:
        for n in self.groups[coords[0]]: 
            for i in self.groups(coords[1]):
                dividend[n+i]=self.differences(n,i)
        return np.mean(dividend)

    def find_min(self):
        '''
        Runs through the table and finds the smallest difference
        Returns a list with the minimum value and coordinates for it in the array [val,coords]
        WARNING: two species with identical sequences will break the process at this step
        '''
        val = np.min(self.grouped_table[np.nonzero(self.grouped_table)])
        return [val,np.where(self.grouped_table == val)]

    def shrink_table(self):
        '''
        Does all necessary steps for the reduction of the table by one cell:
            Finds the minimum value in the table (find_min) (REMINDER, IDENTICAL SEQUENCES WILL BREAK THIS PROCESS)
            Merges the intersecting groups (merge)
            
        '''
        mdata = self.find_min
        min = mdata[0]

        
        return 0
    
    def do_your_thing(self):
        '''
        Executes all methods needed to gather the data needed for the full phylogeny in the following order:
        1. Builds the table (construct_table)
        '''
        self.construct_table()


'''

'''
class PhyTree:

    def __init__(self):
        pass
