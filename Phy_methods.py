import numpy as np
import re
from statistics import mean

def checkSequences(sequences): #Ensures that all AA sequences are valid
    allowed_AAs = ['A','R','N','D','C','Q','E','G','H',
    'I','L','K','M','F','P','O','S','U','T','W','Y','V',]

    for n in range(len(sequences)): #Ensures that all AA sequences are the same length
        if len(sequences[n]) != len(sequences[0]):
            return 'ERROR: ALL STRINGS MUST BE SAME LENGTH'
    
    for n in range(len(sequences)): #Ensures that AA sequences don't contain any incorrect AA 1-letter codes
        for letter in range(len(sequences[0])):
            if sequences[n][letter] not in allowed_AAs:
                return 'ERROR: SEQUENCES(S) CONTAIN INVALID AMINO ACID CODES'

    return "Yep, you're all good :)"

def tableSetup(sequences): #Takes a list of sequences and puts their differences into an array
    
    species_num = len(sequences)
    length = len(sequences[0]) 
    
    comp_table = np.zeros((species_num,species_num)) #sets up an array

    for base in range(species_num):
        for comp in range(species_num):
            count = 0
            for letter in range(len(sequences[base])): #loop finds all different characters 
                if sequences[base][letter] != sequences[comp][letter]:
                    count += 1
            count = round(100*abs(np.log(1-count/length)),2) #accounts for pairwise distance 
            comp_table[base][comp] = count

            for row in range(species_num): #makes tables more readable by removing duplicate data and removing extra 0s
                zero_found = False
                for col in range(species_num):
                    if comp_table[row][col] == 0:
                        zero_found = True
                    if zero_found == True:
                        comp_table[row][col] = 0

    return comp_table

def shrinkTable(givens,groups): #takes values already calculated, and returns the distance between 
    output = []
    givens = np.array(givens)
    groups = list(groups)
    new_row = []
    dimensions = len(givens[0])

    #Determines the smallest difference
    relatives = np.min(givens[np.nonzero(givens)])
    coords = np.where(givens == relatives)

    #assigns indices to the two chosen species
    spec_1_idx = int(re.findall(r'\d+',str(coords))[0])
    spec_2_idx = int(re.findall(r'\d+',str(coords))[1])

    #Extracts the two species' names and gives them a new name as a cluster, while removing them from the original list
    rel_list = [spec_1_idx,spec_2_idx]
    new_name = groups[spec_1_idx] + ', & ' + groups[spec_2_idx]
    for spec in sorted(rel_list, reverse = True):
        del groups[spec]
    groups.append(new_name)

    #Extracts the differences between the two species and all other species
    spec_1r = list(givens[spec_1_idx,:]) #Extracts row 1
    spec_2r = list(givens[spec_2_idx,:]) #Extracts row 2
    spec_1c = list(givens[:,spec_1_idx]) #Extracts column 1
    spec_2c = list(givens[:,spec_2_idx]) #Extracts column 2
    spec_1 = spec_1r + spec_1c
    spec_2 = spec_2r + spec_2c
    
    for n in reversed(range(len(spec_1))):
        if int(spec_1[n]) == 0 or spec_1[n] == relatives :
            spec_1.pop(n)
        if int(spec_2[n]) == 0 or spec_2[n] == relatives :
            spec_2.pop(n)
    
    #Merges the differences, creating a new column

    new_row = []
    for n in range(len(spec_1)):
        pair = [spec_1[n],spec_2[n]]
        average_dif = round(mean(pair),2)
        new_row.append(average_dif)
        


    #Shrinks the table 
    shrink = np.delete(givens,coords[0],0)
    shrink = np.delete(shrink,coords[1],0)
    shrink = np.delete(shrink,coords[0],1)
    shrink = np.delete(shrink,coords[1],1)

    #UNFINISHED: Adds the new row to the bottom and replaces the outermost column

    new_column = []

    for n in range(dimensions-1):
        new_column.append([0])

    shrink = np.row_stack((shrink,[new_row]))
    shrink = np.column_stack((shrink,new_column))


    #Compiles the output
    output.append(shrink)
    output.append(groups)
    output.append(relatives)
    output.append(coords[0])
    output.append(coords[1])

    

    return output

def buildTree(sequences,species):
    species = list(species)
    table = tableSetup(sequences)
    dimensions = table.size
    
    for n in range(dimensions):
        
        shrinkTable()
        


        
    return 0

def generalize(species): #FML

    return 0