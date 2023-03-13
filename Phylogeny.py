import Phy_methods as phy




#Note, AA is short for Amino Acid

phy_file = open('/Users/gabriele/Desktop/Lesiure/Coding/UPGMA_Project/Sequences.txt','r')
phy_string = phy_file.read()
phy_list = phy_string.split('\n')

sequence_key = {} 
sequence_list = []
species_list = []
for n in range(0,len(phy_list)-1,2): #creates a list of the sequences and a dictionary for matching with species names
    sequence_key.update({phy_list[n]:phy_list[n+1]})
    sequence_list.append(phy_list[n+1])
    species_list.append(phy_list[n])

#print(checkSequences(sequence_list))

table = phy.tableSetup(sequence_list)
spec = species_list

print(table)

#table2 = phy.shrinkTable(table,species_list)[0]
#spec2 = phy.shrinkTable(table,species_list)[1]

#print(phy.shrinkTable(table,species_list)[0])

#print(phy.shrinkTable(table2,spec2))


for n in range(4):
    table = phy.shrinkTable(table,spec)[0]
    spec = phy.shrinkTable(table,spec)[1]
    print(table)


phy_file.close()


