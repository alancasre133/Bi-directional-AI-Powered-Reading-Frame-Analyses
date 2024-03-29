import numpy as np 
from collections import defaultdict
from bs4 import BeautifulSoup
import urllib.request
import re
stop_codons = ["tag",'taa','tga']
start_codon = ["atg"]
codons_to_aminoacids = {
  "M" : {
      "nucleotidic_sequence" : ["atg"]
  },
  "P" : {
    "nucleotidic_sequence" : ["ccc","cca","ccg","cct"]
  },
  "A" : {
    "nucleotidic_sequence" : ["gct","gcc","gca","gcg"]
  },
  "R" : {
    "nucleotidic_sequence" : ["cgt","cgc","cga","cgg","aga","agg"]
  },
  "N" : {
    "nucleotidic_sequence" : ["aac","aat"]
  },
  "D" : {
    "nucleotidic_sequence" : ["gat","gac"]
  },
  "C" : {
    "nucleotidic_sequence" : ["tgt","tgc"]
  },
  "E" : {
    "nucleotidic_sequence" : ["gaa","gag"]
  },
  "C" : {
    "nucleotidic_sequence" : ["tgt","tgc"]
  },
  "Q" : {
    "nucleotidic_sequence" : ["caa","cag"]
  },
  "G" : {
    "nucleotidic_sequence" : ["ggt","ggc","gga","ggg"]
  },
  "H" : {
    "nucleotidic_sequence" : ["cat","cac"]
  },
  "I" : {
    "nucleotidic_sequence" : ["att","atc","ata"]
  },
  "L" : {
    "nucleotidic_sequence" : ["tta","ttg","ctt","ctc","cta","ctg"]
  },
  "K" : {
    "nucleotidic_sequence" : ["aaa","aag"]
  },
  "F" : {
    "nucleotidic_sequence" : ["ttt","ttc"]
  },
  "S" : {
    "nucleotidic_sequence" : ["tct","tcc","tca","tcg","agt","agc"]
  },
  "T" : {
    "nucleotidic_sequence" : ["act","acc","aca","acg"]
  },
  "W" : {
    "nucleotidic_sequence" : ["tgg"]
  },
  "Y" : {
    "nucleotidic_sequence" : ["tat","tac"]
  },
  "V" : {
    "nucleotidic_sequence" : ["gtt","gtc","gta","gtg"]
  },
  "Z" : {
    "nucleotidic_sequence" : ["taa","tag","tga"]
  }
}
class NucleotideSequenceAnalyzer():
    def __init__(self) -> None:
        self.codons_to_aminoacids = {
        "M" : {
            "nucleotidic_sequence" : ["atg"]
        },
        "P" : {
            "nucleotidic_sequence" : ["ccc","cca","ccg","cct"]
        },
        "A" : {
            "nucleotidic_sequence" : ["gct","gcc","gca","gcg"]
        },
        "R" : {
            "nucleotidic_sequence" : ["cgt","cgc","cga","cgg","aga","agg"]
        },
        "N" : {
            "nucleotidic_sequence" : ["aac","aat"]
        },
        "D" : {
            "nucleotidic_sequence" : ["gat","gac"]
        },
        "C" : {
            "nucleotidic_sequence" : ["tgt","tgc"]
        },
        "E" : {
            "nucleotidic_sequence" : ["gaa","gag"]
        },
        "C" : {
            "nucleotidic_sequence" : ["tgt","tgc"]
        },
        "Q" : {
            "nucleotidic_sequence" : ["caa","cag"]
        },
        "G" : {
            "nucleotidic_sequence" : ["ggt","ggc","gga","ggg"]
        },
        "H" : {
            "nucleotidic_sequence" : ["cat","cac"]
        },
        "I" : {
            "nucleotidic_sequence" : ["att","atc","ata"]
        },
        "L" : {
            "nucleotidic_sequence" : ["tta","ttg","ctt","ctc","cta","ctg"]
        },
        "K" : {
            "nucleotidic_sequence" : ["aaa","aag"]
        },
        "F" : {
            "nucleotidic_sequence" : ["ttt","ttc"]
        },
        "S" : {
            "nucleotidic_sequence" : ["tct","tcc","tca","tcg","agt","agc"]
        },
        "T" : {
            "nucleotidic_sequence" : ["act","acc","aca","acg"]
        },
        "W" : {
            "nucleotidic_sequence" : ["tgg"]
        },
        "Y" : {
            "nucleotidic_sequence" : ["tat","tac"]
        },
        "V" : {
            "nucleotidic_sequence" : ["gtt","gtc","gta","gtg"]
        },
        "Z" : {
            "nucleotidic_sequence" : ["taa","tag","tga"]
        }
        }
        self.stop_codons = ["tag",'taa','tga']
        self.start_codon = ["atg"]
    
    def __Inrange(aSS,aEE,bSS,bEE):
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
        if (linea_2_inicio >= linea_1_inicio and linea_2_inicio <= linea_1_fin) or (linea_2_fin >= linea_1_inicio and linea_2_fin <= linea_1_fin):
            return 1
        else: 
            return 0


    def translate_nucleotide_sequence_to_aminoacids(self,sequence):
        def convert_nucleotide_sequence_to_aminoacid_sequence(nucleotidic_sequence):
            aminoacid_chain = ""
            for Three in range(0,len(nucleotidic_sequence),3):
                for i in codons_to_aminoacids:
                    if i!='Z' and nucleotidic_sequence[Three:Three+3] in codons_to_aminoacids[i]["nucleotidic_sequence"]:
                        aminoacid_chain += i
                        break
            return aminoacid_chain
        return convert_nucleotide_sequence_to_aminoacid_sequence(sequence)

    def getORFS_and_metadata_fromsequence(self,sequence,minaminoacidlength,include_reverse_orfs=True,include_front_ORFS=True,OrganismName=None,Definition=None):
        def search_proteins_by_index_reverse_circle(full_nucleotide_sequence,start):
            add=0
            reading_frames = []
            reading = "right to left"
            reading_l = []
            start_c=0
            end_c=0
            start_l = []
            end_l = []
            orf_length=0
            ORF_length_l = []
            new_sequence = ""
            read_frame = ""
            for i in full_nucleotide_sequence:
                if i=='a':
                    new_sequence += 't'
                if i=='t':
                    new_sequence += 'a'
                if i=='c':
                    new_sequence += 'g'
                if i=='g':
                    new_sequence += 'c'
            circle_nucleotidic_sequence = 2*''.join(reversed(new_sequence))
            
            for i in range(start, len(circle_nucleotidic_sequence), 3):
                if add==0 and i>(len(circle_nucleotidic_sequence))/2:
                    break
                else:
                    if add==0 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["M"]["nucleotidic_sequence"]:
                        start_c=(len(circle_nucleotidic_sequence)/2)-i
                        add=1
                    if(add==1):
                        read_frame += circle_nucleotidic_sequence[i:i+3]
                    if add==1 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["Z"]["nucleotidic_sequence"]:
                        end_c=(len(circle_nucleotidic_sequence)/2)-i-3
                        start_l.append(start_c)
                        end_l.append(end_c)
                        reading_l.append(reading)
                        reading_frames.append(read_frame)
                        orf_length=start_c-end_c
                        ORF_length_l.append(orf_length)
                        add=0
                        read_frame=""
            return reading_frames,reading_l,start_l,end_l,ORF_length_l
        def search_proteins_by_index_circle(full_nucleotide_sequence,start):
            add=0
            reading_frames = []
            read_frame = ""
            reading = "left to right"
            reading_l = []
            start_c=0
            end_c=0
            start_l = []
            end_l = []
            orf_length=0
            ORF_length_l = []
            circle_nucleotidic_sequence = 2*full_nucleotide_sequence
            for i in range(start, len(circle_nucleotidic_sequence), 3):
                if add==0 and i>(len(circle_nucleotidic_sequence))/2:
                    break
                else:
                    if add==0 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["M"]["nucleotidic_sequence"]:
                        start_c=i
                        add=1
                    if(add==1):
                        read_frame += circle_nucleotidic_sequence[i:i+3]
                    if add==1 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["Z"]["nucleotidic_sequence"]:
                        end_c=i+3
                        start_l.append(start_c)
                        end_l.append(end_c)
                        reading_l.append(reading)
                        reading_frames.append(read_frame)
                        orf_length=end_c-start_c
                        ORF_length_l.append(orf_length)
                        add=0
                        read_frame=""
            return reading_frames,reading_l,start_l,end_l,ORF_length_l
        def getStart_and_stop_codons_front(full_nucleotide_sequence, include_reverse_orfs=True, include_front_ORFS=True):
            add = 0
            reading_frames = []  # Change this to a list
            if include_reverse_orfs == True and include_front_ORFS == True:
                for i in range(0, 3):
                    reading_frames.extend(search_proteins_by_index_reverse_circle(full_nucleotide_sequence, i))
                    reading_frames.extend(search_proteins_by_index_circle(full_nucleotide_sequence, i))
            elif include_front_ORFS == True:
                for i in range(0, 3):
                    reading_frames.extend(search_proteins_by_index_circle(full_nucleotide_sequence, i))
            elif include_reverse_orfs == True:
                for i in range(0, 3):
                    reading_frames.extend(search_proteins_by_index_reverse_circle(full_nucleotide_sequence, i))
            return reading_frames

        list_of_list_ofsequences = getStart_and_stop_codons_front(sequence,include_reverse_orfs,include_front_ORFS)
        metadata = {'ORF':[] , 'Reading Direction': [], 'Start codon': [], 'End codon': [], 'ORFLength':[],'NucleotideSequence':[]}
        for i in range(0,len(list_of_list_ofsequences)//5):
            for x in range(0,len(list_of_list_ofsequences[i*5+0])):
                if len(Transle_nucleotide_sequence_to_aminoacids(list_of_list_ofsequences[i*5+0][x]))>=minaminoacidlength:
                    metadata['ORF'].append(Transle_nucleotide_sequence_to_aminoacids(list_of_list_ofsequences[i*5+0][x]))
                    metadata['Reading Direction'].append(list_of_list_ofsequences[i*5+1][x])
                    metadata['Start codon'].append(int(list_of_list_ofsequences[i*5+2][x]))
                    metadata['End codon'].append(int(list_of_list_ofsequences[i*5+3][x]))
                    metadata['ORFLength'].append(int(list_of_list_ofsequences[i*5+4][x]))
                    metadata['NucleotideSequence'].append(list_of_list_ofsequences[i*5+0][x])
        frames_list=[]
        for x in range(0,len(metadata['End codon'])):
            nesting = 0
            for i in range(0,len(metadata['End codon'])):
                if x!=i and metadata['ORFLength'][x]<metadata['ORFLength'][i] and Inrange(metadata['End codon'][x],metadata['Start codon'][x],metadata['End codon'][i],metadata['Start codon'][i])==1:
                    nesting+=1
            frames_list.append(nesting)
        metadata['Nesting Level'] = frames_list
        Organisms= []
        if OrganismName!=None:
            for nested in metadata["Nesting Level"]:
                Organisms.append(OrganismName)
            metadata['Organism Name']=Organisms
        if Definition!=None:
            Organisms = []
            for nested in metadata["Nesting Level"]:
                Organisms.append(Definition)
            metadata['Definition']=Organisms
        return metadata


def Transle_nucleotide_sequence_to_aminoacids(sequence):
    def convert_nucleotide_sequence_to_aminoacid_sequence(nucleotidic_sequence):
        aminoacid_chain = ""
        for Three in range(0,len(nucleotidic_sequence),3):
            for i in codons_to_aminoacids:
                if i!='Z' and nucleotidic_sequence[Three:Three+3] in codons_to_aminoacids[i]["nucleotidic_sequence"]:
                    aminoacid_chain += i
                    break
        return aminoacid_chain
    return convert_nucleotide_sequence_to_aminoacid_sequence(sequence)
def getORFS_fromsequence(sequence,minaminoacidlength,include_reverse_orfs=True,include_front_ORFS=True):
    def search_proteins_by_index_reverse_circle(full_nucleotide_sequence,start):
        add=0
        reading_frames = []
        new_sequence = ""
        read_frame = ""
        for i in full_nucleotide_sequence:
            if i=='a':
                new_sequence += 't'
            if i=='t':
                new_sequence += 'a'
            if i=='c':
                new_sequence += 'g'
            if i=='g':
                new_sequence += 'c'
        circle_nucleotidic_sequence = 2*''.join(reversed(new_sequence))
        
        for i in range(start, len(circle_nucleotidic_sequence), 3):
            if add==0 and i>(len(circle_nucleotidic_sequence))/2:
                break
            else:
                if add==0 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["M"]["nucleotidic_sequence"]:
                    add=1
                if(add==1):
                    read_frame += circle_nucleotidic_sequence[i:i+3]
                if add==1 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["Z"]["nucleotidic_sequence"]:
                    reading_frames.append(read_frame)
                    add=0
                    read_frame=""
        return reading_frames
    def search_proteins_by_index_circle(full_nucleotide_sequence,start):
        add=0
        reading_frames = []
        read_frame = ""
        circle_nucleotidic_sequence = 2*full_nucleotide_sequence
        for i in range(start, len(circle_nucleotidic_sequence), 3):
            if add==0 and i>(len(circle_nucleotidic_sequence))/2:
                break
            else:
                if add==0 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["M"]["nucleotidic_sequence"]:
                    add=1
                if(add==1):
                    read_frame += circle_nucleotidic_sequence[i:i+3]
                if add==1 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["Z"]["nucleotidic_sequence"]:
                    reading_frames.append(read_frame)
                    add=0
                    read_frame=""
        return reading_frames
    def getStart_and_stop_codons_front(full_nucleotide_sequence,include_reverse_orfs=True,include_front_ORFS=True):
        add = 0
        reading_frames = []
        if include_reverse_orfs==True and include_front_ORFS==True:
            for i in range(0,3):
                reading_frames.append(search_proteins_by_index_reverse_circle(full_nucleotide_sequence,i))
                reading_frames.append(search_proteins_by_index_circle(full_nucleotide_sequence,i))
        elif include_front_ORFS==True:
            for i in range(0,3):
                reading_frames.append(search_proteins_by_index_circle(full_nucleotide_sequence,i))
        elif include_reverse_orfs==True:
            for i in range(0,3):
                reading_frames.append(search_proteins_by_index_reverse_circle(full_nucleotide_sequence,i))
        return reading_frames
    list_of_list_ofsequences = getStart_and_stop_codons_front(sequence,include_reverse_orfs,include_front_ORFS)

    TotalORFS = []
    for list_of_sequences in list_of_list_ofsequences:
        for sequence in list_of_sequences:
            if len(Transle_nucleotide_sequence_to_aminoacids(sequence))>=minaminoacidlength:
                TotalORFS.append(Transle_nucleotide_sequence_to_aminoacids(sequence))  
    return TotalORFS

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
  if (linea_2_inicio >= linea_1_inicio and linea_2_inicio <= linea_1_fin) or (linea_2_fin >= linea_1_inicio and linea_2_fin <= linea_1_fin):
      return 1
  else: 
     return 0


def getORFS_and_metadata_fromsequence(sequence,minaminoacidlength,include_reverse_orfs=True,include_front_ORFS=True,OrganismName=None,Definition=None):
    def search_proteins_by_index_reverse_circle(full_nucleotide_sequence,start):
        add=0
        reading_frames = []
        reading = "right to left"
        reading_l = []
        start_c=0
        end_c=0
        start_l = []
        end_l = []
        orf_length=0
        ORF_length_l = []
        new_sequence = ""
        read_frame = ""
        for i in full_nucleotide_sequence:
            if i=='a':
                new_sequence += 't'
            if i=='t':
                new_sequence += 'a'
            if i=='c':
                new_sequence += 'g'
            if i=='g':
                new_sequence += 'c'
        circle_nucleotidic_sequence = 2*''.join(reversed(new_sequence))
        
        for i in range(start, len(circle_nucleotidic_sequence), 3):
            if add==0 and i>(len(circle_nucleotidic_sequence))/2:
                break
            else:
                if add==0 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["M"]["nucleotidic_sequence"]:
                    start_c=(len(circle_nucleotidic_sequence)/2)-i
                    add=1
                if(add==1):
                    read_frame += circle_nucleotidic_sequence[i:i+3]
                if add==1 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["Z"]["nucleotidic_sequence"]:
                    end_c=(len(circle_nucleotidic_sequence)/2)-i-3
                    start_l.append(start_c)
                    end_l.append(end_c)
                    reading_l.append(reading)
                    reading_frames.append(read_frame)
                    orf_length=start_c-end_c
                    ORF_length_l.append(orf_length)
                    add=0
                    read_frame=""
        return reading_frames,reading_l,start_l,end_l,ORF_length_l
    def search_proteins_by_index_circle(full_nucleotide_sequence,start):
        add=0
        reading_frames = []
        read_frame = ""
        reading = "left to right"
        reading_l = []
        start_c=0
        end_c=0
        start_l = []
        end_l = []
        orf_length=0
        ORF_length_l = []
        circle_nucleotidic_sequence = 2*full_nucleotide_sequence
        for i in range(start, len(circle_nucleotidic_sequence), 3):
            if add==0 and i>(len(circle_nucleotidic_sequence))/2:
                break
            else:
                if add==0 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["M"]["nucleotidic_sequence"]:
                    start_c=i
                    add=1
                if(add==1):
                    read_frame += circle_nucleotidic_sequence[i:i+3]
                if add==1 and circle_nucleotidic_sequence[i:i+3] in codons_to_aminoacids["Z"]["nucleotidic_sequence"]:
                    end_c=i+3
                    start_l.append(start_c)
                    end_l.append(end_c)
                    reading_l.append(reading)
                    reading_frames.append(read_frame)
                    orf_length=end_c-start_c
                    ORF_length_l.append(orf_length)
                    add=0
                    read_frame=""
        return reading_frames,reading_l,start_l,end_l,ORF_length_l
    def getStart_and_stop_codons_front(full_nucleotide_sequence, include_reverse_orfs=True, include_front_ORFS=True):
        add = 0
        reading_frames = []  # Change this to a list
        if include_reverse_orfs == True and include_front_ORFS == True:
            for i in range(0, 3):
                reading_frames.extend(search_proteins_by_index_reverse_circle(full_nucleotide_sequence, i))
                reading_frames.extend(search_proteins_by_index_circle(full_nucleotide_sequence, i))
        elif include_front_ORFS == True:
            for i in range(0, 3):
                reading_frames.extend(search_proteins_by_index_circle(full_nucleotide_sequence, i))
        elif include_reverse_orfs == True:
            for i in range(0, 3):
                reading_frames.extend(search_proteins_by_index_reverse_circle(full_nucleotide_sequence, i))
        return reading_frames

    list_of_list_ofsequences = getStart_and_stop_codons_front(sequence,include_reverse_orfs,include_front_ORFS)
    metadata = {'ORF':[] , 'Reading Direction': [], 'Start codon': [], 'End codon': [], 'ORFLength':[],'NucleotideSequence':[]}
    for i in range(0,len(list_of_list_ofsequences)//5):
        for x in range(0,len(list_of_list_ofsequences[i*5+0])):
            if len(Transle_nucleotide_sequence_to_aminoacids(list_of_list_ofsequences[i*5+0][x]))>=minaminoacidlength:
              metadata['ORF'].append(Transle_nucleotide_sequence_to_aminoacids(list_of_list_ofsequences[i*5+0][x]))
              metadata['Reading Direction'].append(list_of_list_ofsequences[i*5+1][x])
              metadata['Start codon'].append(int(list_of_list_ofsequences[i*5+2][x]))
              metadata['End codon'].append(int(list_of_list_ofsequences[i*5+3][x]))
              metadata['ORFLength'].append(int(list_of_list_ofsequences[i*5+4][x]))
              metadata['NucleotideSequence'].append(list_of_list_ofsequences[i*5+0][x])
    frames_list=[]
    for x in range(0,len(metadata['End codon'])):
        nesting = 0
        for i in range(0,len(metadata['End codon'])):
            if x!=i and metadata['ORFLength'][x]<metadata['ORFLength'][i] and Inrange(metadata['End codon'][x],metadata['Start codon'][x],metadata['End codon'][i],metadata['Start codon'][i])==1:
                nesting+=1
        frames_list.append(nesting)
    metadata['Nesting Level'] = frames_list
    Organisms= []
    if OrganismName!=None:
        for nested in metadata["Nesting Level"]:
            Organisms.append(OrganismName)
        metadata['Organism Name']=Organisms
    if Definition!=None:
        Organisms = []
        for nested in metadata["Nesting Level"]:
            Organisms.append(Definition)
        metadata['Definition']=Organisms
    return metadata

def get_nucleotidesequencefrom_gbfile(filename):
    def getSequence(string):
        pattern = r'ORIGIN\s+([\s\S]*?)//'
        nucleotides_noise = re.findall(pattern,string,re.DOTALL)
        pattern = r'[acgt]+'
        full_nucleotide_sequence = ""
        full_nucleotide_sequences = []
        for every in nucleotides_noise:
            nucleotides = re.findall(pattern,every,re.DOTALL)
            for nucleotide_sequence in nucleotides:
                full_nucleotide_sequence += nucleotide_sequence
            full_nucleotide_sequences.append(full_nucleotide_sequence)
            full_nucleotide_sequence = ""
        return full_nucleotide_sequences

    #Obtiene las secuencias de nucleotidos completas por cada virus en el archivo
    path = r"assets\\"
    file = open(path+filename,"r")
    Data = file.read()
    Data = getSequence(Data)
    return Data

def get_nucleotidesequence_and_metadata_from_gbfile(filename):
    def getSequence(string):
        pattern = r'ORIGIN\s+([\s\S]*?)//'
        nucleotides_noise = re.findall(pattern,string,re.DOTALL)
        pattern = r'[acgt]+'
        full_nucleotide_sequence = ""
        full_nucleotide_sequences = []
        for every in nucleotides_noise:
            nucleotides = re.findall(pattern,every,re.DOTALL)
            for nucleotide_sequence in nucleotides:
                full_nucleotide_sequence += nucleotide_sequence
            full_nucleotide_sequences.append(full_nucleotide_sequence)
            full_nucleotide_sequence = ""
        return full_nucleotide_sequences
    def getOrganism(string):
        pattern = r'ORGANISM\s+([\s\S]*?)\n'
        nucleotides_noise = re.findall(pattern,string,re.DOTALL)
        return nucleotides_noise
    def getDefinition(string):
        pattern = r'DEFINITION\s+([\s\S]*?)ACCESSION'
        nucleotides_noise = re.findall(pattern,string,re.DOTALL)
        new = []
        a = ""
        for i in nucleotides_noise:
            a= i.replace("\n", "")
            a = a.replace("  ","")
            new.append(a.strip())
        return new
    #Obtiene las secuencias de nucleotidos completas por cada virus en el archivo
    path = r"assets\\"
    file = open(path+filename,"r")
    Data1 = file.read()
    Data = getSequence(Data1)

    return Data,getOrganism(Data1),getDefinition(Data1)


def returninsertion(organism1,organism2,startindex=3):
    dictio = defaultdict(int)
    res=""
    for i in organism2:
        dictio[i]+=1
    stop = defaultdict(int)
    for i in organism2[len(organism2)-3:len(organism2)]:
        dictio[i]-=1
        stop[i]+=1
    for i in organism1:
        n = ((dictio[i]-1)>=0)
        if(n):
            dictio[i]-=1
            res+=i
        else:
            if((dictio[i]-1)<0):
                return False
    val=[]
    for v in dictio.keys():
        val.append([v,dictio[v]])
    print(val)
    print(res)
    while val:
        idx=np.random.randint(low=0,high=len(val))
        if val[idx][1]>0:
            res+=val[idx][0]
            val[idx][1]-=1
        else:
            val[idx]=val[-1]
            val.pop()
    val=[]
    print(stop)
    for v in stop.keys():
        res+=v*stop[v]
    return res
print(returninsertion("atgcfatgatga","atgtgacfatgatga",0))