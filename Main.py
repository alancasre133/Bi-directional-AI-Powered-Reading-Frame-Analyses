from align import *
from BioInformatic import *
from BioInformaticQuery import *
import pandas as pd

geneticcode,Organism,Definition=get_nucleotidesequence_and_metadata_from_gbfile("PorcineCircovirus3.gb")
geneticcode=get_nucleotidesequencefrom_gbfile("PorcineCircovirus3.gb")
print(geneticcode)
index = 0
for j in geneticcode:
  ORFSL = getORFS_and_metadata_fromsequence(j,minaminoacidlength=25,OrganismName=Organism[index],Definition=Definition[index])
  Data = pd.DataFrame(ORFSL)
  Data = Data.reindex(columns=['Definition','Organism Name','NucleotideSequence','ORF','Reading Direction','Nesting Level','Start codon','End codon','ORFLength'])
  Data['full_sequence']=j
  addOrganism_to_database(Data,"DatabaseORFS\\out.csv")
  index+=1
analyzeORF()
#dataset=pd.read_excel("AnalyzedORFS\exported_json_data.xlsx")
