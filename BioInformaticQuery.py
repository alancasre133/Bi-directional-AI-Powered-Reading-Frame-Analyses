import pandas as pd
import os
import csv
import time
import platform

def delimit(ORF_ROWS):
  available =  ORF_ROWS[((ORF_ROWS['Reading Direction'] == 'left to right') & (ORF_ROWS['Nesting Level'] == '0'))]
  Data = ORF_ROWS[((ORF_ROWS['Reading Direction'] == 'left to right') & (ORF_ROWS['Nesting Level'] == '0')) | ((ORF_ROWS['Nesting Level'] != '0') & (ORF_ROWS['Reading Direction'] == 'right to left'))]
  print("These are the ORFS you can analyse, select one to choose his nested ORF to analyse:")
  index=1
  for i in list(available['ORF']):
    print(str(index) + ". " + "Organism Definition: "+ list(available["Definition"])[0] + " , ORF:"+ i)
    index+=1
def analyzeoption(rows,startquery,endquery,columnnames):
  print()  
  Data = pd.DataFrame(rows[startquery:endquery],columns=columnnames)
  print("Advertisement: For the analysis of ORFS in this algorithm we will only include those that go in the direction of reading from left to right and the ORFS nested to this")
  print(Data)
  delimit(Data)
def cleanscreen():
  time.sleep(1)
  if platform.system()=='Windows':
    os.system('cls')
  else:
    os.system('clear')

def OrganismInDataBase(path,Organism):
  with open(path, "r") as entrada:
    csv_reader = csv.reader(entrada, delimiter=',')
    for fila in csv_reader:
      if fila[0] == Organism:
        return 1
  return 0
def addOrganism_to_database(OrganismMetadata,path):
  OrganismMetadata = OrganismMetadata
  if (os.path.isfile(path) and os.stat(path).st_size != 0):
    if OrganismInDataBase(path,OrganismMetadata['Definition'][0])!=1:
      OrganismMetadata.to_csv("DatabaseORFS\\out.csv",index=False,mode="a",header=False)
  elif (os.path.isfile(path)):
    OrganismMetadata.to_csv("DatabaseORFS\\out.csv",index=False,mode="w")
  else:
    os.makedirs(path, exist_ok=True)  
def reviewgeneticodes(path="DatabaseORFS\\out.csv"):
  Organismslist = []
  OrfsList = []
  orfcount = 0
  first =True
  row_list = []
  with open(path, "r") as entrada:
    csv_reader = csv.reader(entrada, delimiter=',')
    for fila in csv_reader:
      orfcount+=1
      row_list.append(fila)
      if fila[0] not in Organismslist and first==False:
        Organismslist.append(fila[0])
        if len(Organismslist)==1:
          orfcount=1
        else:
          OrfsList.append(orfcount-1)
          orfcount=1
      else:
        first=False
        
    OrfsList.append(orfcount)
  index=1
  print("Choose a Organism to analyze their ORF:")
  for Organism in Organismslist:
    print(str(index)+ ". " + Organism + " ; ORFS: ",str(OrfsList[index-1]))
    index+=1
  op=int(input())
  while op>index-1 or op<1:
    print("Thats not a correct option, please try again")
    op=int(input())
  startquery=1
  endquery=0
  for i in range(0,len(Organismslist)):
    if i==op-1:
      endquery=startquery+OrfsList[i]
      break
    startquery+=OrfsList[i]    
  cleanscreen()
  analyzeoption(row_list,startquery=startquery,endquery=endquery,columnnames=row_list[0])

def analyzeORF():
  print("Choose an option:")
  print("1. Review ORFS being analyzed")
  print("2. Review Genetic Codes and analyze a specific ORF")
  cleanscreen()
  reviewgeneticodes()