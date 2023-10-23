import pandas as pd
import os

def addOrganism_to_database(OrganismMetadata,path):
  OrganismMetadata = OrganismMetadata
  if (os.path.isfile(path) and os.stat(path).st_size != 0):
    OrganismMetadata.to_csv("DatabaseORFS\\out.csv",index=False,mode="a",header=False)
  elif (os.path.isfile(path)):
    OrganismMetadata.to_csv("DatabaseORFS\\out.csv",index=False,mode="w")
  else:
    os.makedirs(path, exist_ok=True)  