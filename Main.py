from align import *
from BioInformatic import *
from BioInformaticQuery import *
import pandas as pd
geneticcode,Organism=get_nucleotidesequence_and_metadata_from_gbfile("PorcineCircovirus3.gb")
geneticcode=get_nucleotidesequencefrom_gbfile("PorcineCircovirus3.gb")
for j in geneticcode:
  ORFSL = getORFS_and_metadata_fromsequence(j,minaminoacidlength=25,OrganismName=Organism)
  Data = pd.DataFrame(ORFSL)
  Data = Data.reindex(columns=['Organism Name','NucleotideSequence','ORF','Reading Direction','Nesting Level','Start codon','End codon','ORFLength'])
  addOrganism_to_database(Data,"DatabaseORFS\\out.csv")