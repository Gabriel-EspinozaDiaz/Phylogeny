import numpy as np
import math
import re

'''
Reads in data and constructs a pairwise matrix
Subsequently constructs a table by Constructs a table using unweighted pairwise distance by methods of arithmetic means (UPGMA)
- 
'''
class Pairwise:
    def __init__(self,directory):
        if type(directory) != str:
            raise TypeError('DIRECTORY MUST BE IN STRING FORM')
        phy_file = open(directory,'r')
        phy_string = phy_file.read()
        phy_list = phy_string.split('\n')
        phy_file.close()

        #Generates a dictionary with species names as keys and sequences as items (specDict), 
        #with a parent list for referencing (specList), ensuring that all inputs are correct.
        #Additionally stores a copy of the list in self.groups that will store species as clades during UPGMA
        self.spec_dict = {} 
        self.specList = [] 
        for n in range(0,len(phy_list),2): 
            self.specList.append(phy_list[n])
            self.spec_dict.update({phy_list[n]:phy_list[n+1]})
        self.specNum = len(self.spec_dict) 
        self.check_sequences()
        self.groups = self.specList.copy() 

        #Builds the starting table, comparing the lengths of different species
        self.parent_table = np.zeros((self.specNum,self.specNum),dtype=float) 
        exc = 0
        for n in range(len(self.specList)):
            exc += 1
            for i in range(exc,len(self.specList)):
                self.parent_table[n][i] = self.differences(n,i)
        self.grouped_table = self.parent_table.copy()

    def check_sequences(self):
        '''
        Ensures that the self.specDict only contains valid amino acid sequences. All must be the same length
        '''
        sequences = []
        for n in range(self.specNum):
            sequences.append(self.spec_dict[self.specList[n]])
        allowed_AAs = ['A','R','N','D','C','Q','E','G','H',
    'I','L','K','M','F','P','O','S','U','T','W','Y','V',]
        
        lenSeqs = []
        for n in range(len(sequences)): #Ensures that all AA sequences are the same length
            if len(sequences[n]) != len(sequences[0]):
                    lenSeqs.append(sequences[n]) 
    
        codeSeqs = []
        for n in range(len(sequences)): #Ensures that AA sequences don't contain any incorrect AA 1-letter codes
            for letter in range(len(sequences[0])):
                if sequences[n][letter] not in allowed_AAs:
                    codeSeqs.append(sequences[n])
                    break

        if len(lenSeqs) != 0:
            raise TypeError(f'ALL SEQUENCES MUST BE SAME LENGTH, THE FOLLOWING SEQUENCES HAVE A DIFFERENT NUMBER AMINO ACIDS: {lenSeqs}')
        elif len(codeSeqs) != 0:
            raise TypeError(f'SEQUENCES(S): {codeSeqs} CONTAIN INVALID AMINO ACID CODES')
        else:
            print('ALL VALID AMINO ACID SEQUENCES')

    def differences(self,spec1,spec2):
        '''
        Takes the index of 2 species and calculates and returns the number of pairwise differences.
        Takes multiple hits into account using the Poisson Formula
        '''
        spec1 = self.spec_dict[self.specList[spec1]]
        spec2 = self.spec_dict[self.specList[spec2]]
        count = 0
        for n in range(len(spec1)):
            if spec1[n] != spec2[n]:
                count += 1
        return round(math.log(1-(count/len(spec1)))*-100,2)
        
    def find_min(self):
        '''
        Runs through the table and finds the smallest difference
        Returns a list with the minimum value and coordinates for it in the array [val,coords]
        WARNING: two species with identical sequences will break the process at this step
        '''
        val = np.min(self.grouped_table[np.nonzero(self.grouped_table)])
        return [val,np.where(self.grouped_table == val)]

    def collect_children(self,node):
        children = []
        visits = [node]
        while len(visits) > 0:
            children.append(node.children)
            node = 'idk'
        return children

    def merge(self,coords):
        '''
        Uses the parent table and index to average 2 species/group differences
        Accepts 2 indices for groups on the current table
        Returns the merged value
        '''
        gs1 = coords[0]
        gs2 = coords[1]
        dividend = np.empty(len(gs1)*len(gs2))
        for n in self.groups[coords[0]]: 
            for i in self.groups(coords[1]):
                dividend[n+i]=self.differences(n,i)
        return np.mean(dividend)

    def shrink_table(self):
        '''
        Does all necessary steps for the reduction of the table by one cell:
            Finds the minimum value in the table (find_min) (REMINDER, IDENTICAL SEQUENCES WILL BREAK THIS PROCESS)
            Merges the intersecting groups (merge)
            
        '''
        mdata = self.find_min
        min = mdata[0]

        return 0


'''
Stores a name, a distance 
'''
class PhyNode:
    def __init__(self,name,distance):
        self.name = name
        self.dist = distance
        self.children = []
    
    def add_child(self,child):
        self.children.append(child)
    
    def get_children(self):
        visits = [self.name] + self.children
        node_children = []
        while len(visits) != 0:
            current_node = visits.pop()
            node_children.append(current_node.name)
        return node_children
    