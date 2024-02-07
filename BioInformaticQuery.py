import pandas as pd
import os
import csv
import time
import platform
import json as j
from align import *
class tree:
  def __init__(self) -> None:
    pass
"""
Check if:
   1. A path doesnt exist, file doesnt exist, create the path and the file specified in the path
   2. A path exist, the file doesnt exist, return 0
   3. A path exist, file exist and is empty (doesnt have content), return 1
   4. A path exist, file exist and have content, return 2
"""
def createdjson(path):
  path = path.split("\\")
  value = -1
  if (os.path.isfile(path[0]+"\\"+path[1]) and os.stat(path[0]+"\\"+path[1]).st_size != 0):
    value=2
  elif (os.path.isdir(path[0]) and os.path.isfile(path[0]+"\\"+path[1]) and os.stat(path[0]+"\\"+path[1]).st_size == 0):
    value=1
  elif (os.path.isdir(path[0]) and not os.path.isfile(path[0]+"\\"+path[1])):
    value=0
  else:
    os.makedirs(path[0], exist_ok=True)
  return value
"""

"""
def updatepositions(aSS,aEE,bSS,bEE):
  if aSS>aEE:
     aS=aEE
  elif aSS<=aEE:
     aS=aSS
  if bSS>bEE:
     bS=bEE
     bE=bSS
  elif bSS<=bEE:
     bS=bSS
     bE=bEE
  return bE-aS,bS-aS
"""
Check where a interval is inside of another, like [1,5] [2,3] if one is inside of the other return 1 else return 0
"""
def Inrange(aSS,aEE,bSS,bEE):
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
  if linea_1_inicio <= linea_2_inicio and linea_2_fin <= linea_1_fin:
      return 1
  elif (linea_2_inicio >= linea_1_inicio and linea_2_inicio <= linea_1_fin) or (linea_2_fin >= linea_1_inicio and linea_2_fin <= linea_1_fin):
      return 1
  else: 
     return 0
"""

"""
def createJSON(Data_rows,op):
  print(Data_rows,op)
  if createdjson("AnalyzedORFS\\matchs.json")==2 and createdjson("AnalyzedORFS\\exported_json_data.xlsx")!=1:
    with open('AnalyzedORFS\\matchs.json','r') as json_file:
      d = j.load(json_file)
    df = pd.DataFrame(d)
    excel =pd.read_excel("AnalyzedORFS\\exported_json_data.xlsx")
    index=0
    names=['SpecieName','Ancestor','ORFSequenceNucleotide','ORFSequenceAminoacid','Direction','Startcodon','Endcodon','Startaminoacid','Endaminoacid','AncestorSequence']                                                                                                                                                                                 
    Row_Selected = (Data_rows.iloc[op-1]).to_dict()
    Row_Selected['Similarity']=[]
    Row_Selected['SpecieName']=[]
    Row_Selected['Ancestor']=[]
    Row_Selected['ORFSequenceNucleotide']=[]
    Row_Selected['ORFSequenceAminoacid']=[]
    Row_Selected['Startcodon']=[]
    Row_Selected['Endcodon']=[]
    Row_Selected['Startaminoacid']=[]
    Row_Selected['Endaminoacid']=[]
    Row_Selected['AncestorSequence']=[]
    Row_Selected['Direction']=[]
    for i in range(0,len(Data_rows)):
      if i!=(op-1):
        if len(str(Data_rows.iloc[i]['ORF']))<len(str(Data_rows.iloc[op-1]['ORF'])):
          x=pair_alignment(str(Data_rows.iloc[i]['ORF']),str(Data_rows.iloc[op-1]['ORF']))
        else:
          x=pair_alignment(str(Data_rows.iloc[op-1]['ORF']),str(Data_rows.iloc[i]['ORF']))
        x=x/len(str(Data_rows.iloc[op-1]['ORF']))
        if x>=0.35:
          index=0
          print(f"Match Percentage: {x*100}")
          for m in names:
            Row_Selected[m].append(list(Data_rows.iloc[i])[index])
            index+=1
          Row_Selected['Similarity'].append(str(x))
    with open('DatabaseORFS\\out.csv') as file:
      DataBaseORFS = [fila for fila in csv.reader(file, delimiter=',')]
    os.system('cls')
    Data = pd.DataFrame(DataBaseORFS[1:],columns=DataBaseORFS[0])
    Data_rowssave=Data_rows
    Data_rows=Data
    for i in range(0,len(Data)):
      if i!=(op-1):
        if len(str(Data_rows.iloc[i]['ORF']))<len(str(Data_rowssave.iloc[op-1]['ORF'])):
          x=pair_alignment(str(Data_rows.iloc[i]['ORF']),str(Data_rowssave.iloc[op-1]['ORF']))
        else:
          x=pair_alignment(str(Data_rowssave.iloc[op-1]['ORF']),str(Data_rows.iloc[i]['ORF']))
        x=x/len(str(Data_rowssave.iloc[op-1]['ORF']))
        if x>=0.35:
          print(f"Match Percentage: {x*100}")
          index=0
          for m in names:
            Row_Selected[m].append(list(Data_rows.iloc[i])[index])
            index+=1
          Row_Selected['Similarity'].append(str(x))
    with open('AnalyzedORFS\\matchs.json', 'w') as json_file:
      j.dump(Row_Selected, json_file)
    if createdjson("AnalyzedORFS\\matchs.json")==2 and createdjson("AnalyzedORFS\\exported_json_data.xlsx")==2: 
      with open('AnalyzedORFS\\matchs.json','r') as json_file:
        d = j.load(json_file)
      df = pd.DataFrame(d)
    df_combinado = pd.concat([excel, df], ignore_index=False)
    df_combinado = df_combinado.loc[:, ~df_combinado.columns.str.contains('^Unnamed')]
    df_combinado.to_excel('AnalyzedORFS\\exported_json_data.xlsx')
    return None
  if createdjson("AnalyzedORFS\\matchs.json")==2 and createdjson("AnalyzedORFS\\exported_json_data.xlsx")==1:
    with open('AnalyzedORFS\\matchs.json','r') as json_file:
      d = j.load(json_file)
    df=pd.DataFrame(d)
    df.to_excel('AnalyzedORFS\\exported_json_data.xlsx')
  if createdjson("AnalyzedORFS\\matchs.json")==1:
    names=['SpecieName','Ancestor','ORFSequenceNucleotide','ORFSequenceAminoacid','Direction','Startcodon','Endcodon','Startaminoacid','Endaminoacid','AncestorSequence']                                              
    Row_Selected = (Data_rows.iloc[op-1]).to_dict()
    Row_Selected['Similarity']=[]
    Row_Selected['SpecieName']=[]
    Row_Selected['Ancestor']=[]
    Row_Selected['ORFSequenceNucleotide']=[]
    Row_Selected['ORFSequenceAminoacid']=[]
    Row_Selected['Startcodon']=[]
    Row_Selected['Endcodon']=[]
    Row_Selected['Startaminoacid']=[]
    Row_Selected['Endaminoacid']=[]
    Row_Selected['AncestorSequence']=[]
    Row_Selected['Direction']=[]
    for i in range(0,len(Data_rows)):
      if i!=(op-1):
        if len(str(Data_rows.iloc[i]['ORF']))<len(str(Data_rows.iloc[op-1]['ORF'])):
          x=pair_alignment(str(Data_rows.iloc[i]['ORF']),str(Data_rows.iloc[op-1]['ORF']))
        else:
          x=pair_alignment(str(Data_rows.iloc[op-1]['ORF']),str(Data_rows.iloc[i]['ORF']))
        x=x/len(str(Data_rows.iloc[op-1]['ORF']))
        if x>=0.35:
          print(f"Match Percentage: {x*100}")
          index=0
          for m in names:
            Row_Selected[m].append(list(Data_rows.iloc[i])[index])
            index+=1
          Row_Selected['Similarity'].append(str(x))
    with open('DatabaseORFS\\out.csv') as file:
      DataBaseORFS = [fila for fila in csv.reader(file, delimiter=',')]
    os.system('cls')
    Data = pd.DataFrame(DataBaseORFS[1:],columns=DataBaseORFS[0])
    Data_rowssave=Data_rows
    Data_rows=Data
    for i in range(0,len(Data)):
      if i!=(op-1):
        if len(str(Data_rows.iloc[i]['ORF']))<len(str(Data_rowssave.iloc[op-1]['ORF'])):
          x=pair_alignment(str(Data_rows.iloc[i]['ORF']),str(Data_rowssave.iloc[op-1]['ORF']))
        else:
          x=pair_alignment(str(Data_rowssave.iloc[op-1]['ORF']),str(Data_rows.iloc[i]['ORF']))
        x=x/len(str(Data_rowssave.iloc[op-1]['ORF']))
        if x>=0.35:
          print(f"Match Percentage: {x*100}")
          index=0
          for m in names:
            Row_Selected[m].append(list(Data_rows.iloc[i])[index])
            index+=1
          Row_Selected['Similarity'].append(str(x))
    with open('AnalyzedORFS\\matchs.json', 'w') as json_file:
      j.dump(Row_Selected, json_file)
def getNameFile(Organism,ORFaminoacid):
  OrganismFile = Organism.replace(" ","_")
  Prefix = ORFaminoacid[:3]
  Sufix = ORFaminoacid[(len(ORFaminoacid))//2:(len(ORFaminoacid)//2+3)]
  Pofix = ORFaminoacid[(len(ORFaminoacid))-3:(len(ORFaminoacid))]
  return (OrganismFile+"-"+Prefix+Sufix+Pofix)
"""
From the ORFS in out.csv, will be selected the ones that fulfill the next conditions:
  1. Is nested or is 'Nested' is greater than 0
  2. Direction is right to left
"""
def delimit2(ORF_ROWS):
  available =  ORF_ROWS[((ORF_ROWS['Reading Direction'] == 'left to right') & (ORF_ROWS['Nesting Level'] == '0'))]
  Data = ORF_ROWS[((ORF_ROWS['Reading Direction'] == 'left to right') & (ORF_ROWS['Nesting Level'] == '0')) | ((ORF_ROWS['Nesting Level'] != '0') & (ORF_ROWS['Reading Direction'] == 'right to left'))]
  print("These are the ORFS you can analyse, select one to choose his nested ORF to analyse:")
  
  index=0
  ava=Data.copy()
  for j in list(available['ORF']):
    
    Data=ava.copy()
    print(str(index+1) + ". " + "Organism Definition: "+ list(available["Definition"])[0] + " , ORF:"+ j)
    index+=1 
    op=index
    print(list(available['ORF']), op-1)
    delete=list(available['ORF'])[op-1] 
    kernel = Data['ORF'] == delete
    df = Data.drop(Data[kernel].index)
    kernel = []
    for i in range(0,len(df)):
      kernel.append(not Inrange(int((df.iloc[i])['End codon']),int((df.iloc[i])['Start codon']),int(list(available['End codon'])[op-1]),int(list(available['Start codon'])[op-1]))==1)
    df1 = (df.drop(df[kernel].index))
    for inde,row in df1.iterrows():
      df1.at[inde,'Start codon'],df1.at[inde,'End codon']=updatepositions(int(list(available['End codon'])[op-1]),int(list(available['Start codon'])[op-1]),int(df1.at[inde,'End codon']),int(df1.at[inde,'Start codon']))
    ind=1
    h=0
    cleanscreen()
    for i in range(0,len(df1)):
      print(str(ind) +". "+ str((df1.iloc[i])['Definition']) + " ; ORF sequence: " + str((df1.iloc[i])['ORF']) + " ; OrfLength nucleotides : "   + str((df1.iloc[i])['ORFLength']) + " ; Start Codon : " + str((df1.iloc[i])['Start codon']) + " ; End Codon : " +str((df1.iloc[i])['End codon']) )
      ind+=1
      print("choose one to gather data and analyze among other ORFS in genetic codes")
      opp = i
      createJSON(df1,opp)

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
  delete=list(available['ORF'])[op-1] 
  kernel = Data['ORF'] == delete
  df = Data.drop(Data[kernel].index)
  kernel = []
  for i in range(0,len(df)):
    kernel.append(not Inrange(int((df.iloc[i])['End codon']),int((df.iloc[i])['Start codon']),int(list(available['End codon'])[op-1]),int(list(available['Start codon'])[op-1]))==1)
  df1 = (df.drop(df[kernel].index))
  for index,row in df1.iterrows():
    df1.at[index,'Start codon'],df1.at[index,'End codon']=updatepositions(int(list(available['End codon'])[op-1]),int(list(available['Start codon'])[op-1]),int(df1.at[index,'End codon']),int(df1.at[index,'Start codon']))
  index=1
  h=0
  while True:
    cleanscreen()
    for i in range(0,len(df1)):
      print(str(index) +". "+ str((df1.iloc[i])['Definition']) + " ; ORF sequence: " + str((df1.iloc[i])['ORF']) + " ; OrfLength nucleotides : "   + str((df1.iloc[i])['ORFLength']) + " ; Start Codon : " + str((df1.iloc[i])['Start codon']) + " ; End Codon : " +str((df1.iloc[i])['End codon']) )
      index+=1
    print(str(index) + ". Return")
    print("choose one to gather data and analyze among other ORFS in genetic codes")
    op = int(input())
    if(op-1>=0 and op-1<index):
      cleanscreen()
      h=op
      break
    else:
      print("Choose a valid option")
  if (h!=index):
    print(getNameFile(str((df1.iloc[op-1])['Definition']),str((df1.iloc[op-1])['ORF'])))
    createJSON(df1,op)
def analyzeoption(rows,startquery,endquery,columnnames):
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
  analyze2(row_list,startquery=startquery,endquery=endquery,columnnames=row_list[0])

def analyze2(rows,startquery,endquery,columnnames):
  Data = pd.DataFrame(rows[startquery:endquery],columns=columnnames)
  print("Advertisement: For the analysis of ORFS in this algorithm we will only include those that go in the direction of reading from left to right and the ORFS nested to this")
  print(Data)
  delimit2(Data)

def thirdoption(path="DatabaseORFS\\out.csv"):
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
  print("Organisms to analyze their ORF:")
  for Organism in Organismslist:
    print(str(index)+ ". " + Organism + " ; ORFS: ",str(OrfsList[index-1]))
    index+=1
  startquery=1
  endquery=0
  print(row_list)
  
  for i in range(0,len(Organismslist)):  
    endquery=startquery+OrfsList[i]
    cleanscreen()
    analyze2(row_list,startquery=startquery,endquery=endquery,columnnames=row_list[0])
    startquery+=OrfsList[i]

def analyzeORF():
  print("Choose an option:")
  print("1. Review ORFS being analyzed")
  print("2. Review Genetic Codes and analyze a specific ORF")
  print("3. Analyze all in the Database")
  cleanscreen()
  thirdoption()
  #reviewgeneticodes()
