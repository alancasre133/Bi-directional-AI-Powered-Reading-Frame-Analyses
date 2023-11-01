import pandas as pd
import os
import csv
import time
import platform
import json as j

def createjson(path):
  path = path.split("\\")
  value = -1
  if (os.path.isfile(path[0]+"\\"+path[1]) and os.stat(path[0]+"\\"+path[1]).st_size != 0):
    value=2
  elif (os.path.isdir(path[0]) and os.path.isfile(path[0]+"\\"+path[1]) and os.stat(path[0]+"\\"+path[1]).st_size == 0):
    value=1
  elif (os.path.isdir(path[0]) and not os.path.isfile(path[0]+"\\"+path[1])):
    value=0
    with open(path[0]+"\\"+path[1], 'w') as f:
      print()
  else:
    os.makedirs(path[0], exist_ok=True)
  return value

def updatepositions(aSS,aEE,bSS,bEE):
  if aSS>aEE:
     aS=aEE
     aE=aSS
  elif aSS<=aEE:
     aS=aSS
     aE=aEE
  if bSS>bEE:
     bS=bEE
     bE=bSS
  elif bSS<=bEE:
     bS=bSS
     bE=bEE
  linea_1_inicio = bS
  linea_1_fin = bE

  linea_2_inicio = aS
  linea_2_fin = aE
  return bE-aS,bS-aS
  # Comprobación de si la segunda línea está dentro del rango de la primera línea
  
def Inrange(aSS,aEE,bSS,bEE):
  # Definición de los intervalos
  if aSS>aEE:
     aS=aEE
     aE=aSS
  elif aSS<=aEE:
     aS=aSS
     aE=aEE
  if bSS>bEE:
     bS=bEE
     bE=bSS
  elif bSS<=bEE:
     bS=bSS
     bE=bEE
  linea_1_inicio = bS
  linea_1_fin = bE

  linea_2_inicio = aS
  linea_2_fin = aE
  # Comprobación de si la segunda línea está dentro del rango de la primera línea
  if linea_1_inicio <= linea_2_inicio and linea_2_fin <= linea_1_fin:
      return 1
  elif (linea_2_inicio >= linea_1_inicio and linea_2_inicio <= linea_1_fin) or (linea_2_fin >= linea_1_inicio and linea_2_fin <= linea_1_fin):
      return 1
  else: 
     return 0


def delimit(ORF_ROWS):
  available =  ORF_ROWS[((ORF_ROWS['Reading Direction'] == 'left to right') & (ORF_ROWS['Nesting Level'] == '0'))]
  Data = ORF_ROWS[((ORF_ROWS['Reading Direction'] == 'left to right') & (ORF_ROWS['Nesting Level'] == '0')) | ((ORF_ROWS['Nesting Level'] != '0') & (ORF_ROWS['Reading Direction'] == 'right to left'))]
  print("These are the ORFS you can analyse, select one to choose his nested ORF to analyse:")
  while True:
    index=1
    
    for i in list(available['ORF']):
      print(str(index) + ". " + "Organism Definition: "+ list(available["Definition"])[0] + " , ORF:"+ i)
      index+=1 
    op=int(input())
    if op<index and op>0:
      break
    else:
      print("Choose a correct option")
  #print(list(available['ORF'])[op-1] + " : you choose this")  
  delete=list(available['ORF'])[op-1] 
  #DATA2 = Data.drop((list(available['ORF'])[op-1]))
  
  #print(len(DATA2))
  kernel = Data['ORF'] == delete
  df = Data.drop(Data[kernel].index)
  kernel = []
  for i in range(0,len(df)):
    #Inrange((available['End codon'])[op-1],(available['Start codon'])[op-1],(df.iloc[i])['End codon'][i],(df.iloc[i])['Start codon'][i])==1
    kernel.append(not Inrange(int((df.iloc[i])['End codon']),int((df.iloc[i])['Start codon']),int(list(available['End codon'])[op-1]),int(list(available['Start codon'])[op-1]))==1)
  df1 = (df.drop(df[kernel].index))
  

    #print(str(int(list(available['End codon'])[op-1])) + "," + str(int(list(available['Start codon'])[op-1])) + "," + str(int((df.iloc[i])['End codon']))+ "," + str(int((df.iloc[i])['Start codon'])))
  for index,row in df1.iterrows():
    df1.at[index,'Start codon'],df1.at[index,'End codon']=updatepositions(int(list(available['End codon'])[op-1]),int(list(available['Start codon'])[op-1]),int(df1.at[index,'End codon']),int(df1.at[index,'Start codon']))
  index=1
  while True:
    cleanscreen()
    for i in range(0,len(df1)):
      print(str(index) +". "+ str((df1.iloc[i])['Definition']) + " ; ORF sequence: " + str((df1.iloc[i])['ORF']) + " ; OrfLength nucleotides : "   + str((df1.iloc[i])['ORFLength']) + " ; Start Codon : " + str((df1.iloc[i])['Start codon']) + " ; End Codon : " +str((df1.iloc[i])['End codon']) )
      index+=1
    print(str(index) + ". Return")
    print("choose one to gather data and analyze among other ORFS in genetic codes")
    op = int(input())
    if(op-1>=0 and op-1<index):
      break
    else:
      print("Choose a valid option")
    cleanscreen()
    #imprimir a su padre

# Usar el método drop para eliminar las filas que cumplan con la máscar
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