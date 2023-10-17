import random
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import pandas as pd
import re
import fileinput
import os
import math
import itertools as it
from posixpath import split 
from bs4 import BeautifulSoup
import urllib.request
import re
import sys
from align import *
from BioInformatic import *

def convert_nucleotide_sequence_to_aminoacid_sequence(nucleotidic_sequence):
  aminoacid_chain = ""
  for Three in range(0,len(nucleotidic_sequence),3):
    for i in codons_to_aminoacidsconversion():
      if i!='Z' and nucleotidic_sequence[Three:Three+3] in codons_to_aminoacidsconversion()[i]["nucleotidic_sequence"]:
        aminoacid_chain += i
        break
  return aminoacid_chain

def convert_aminoacid_sequence_to_nucleotide_sequence(aminoacid_sequence):
  nucleotide_chain = ""
  for Three in aminoacid_sequence:
    nucleotide_chain += codons_to_aminoacidsconversion()[Three]["nucleotidic_sequence"][0]
    
  return nucleotide_chain

def limits():
    sys.setrecursionlimit(2000)
    sys.set_int_max_str_digits(9000)   
#A283051
def codons_to_aminoacidsconversion():  
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
    return codons_to_aminoacids

def getSequencefromfile(string):
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

def getDatafromfile(filenamegb):
    filename= filenamegb
    file = open(filename,"r")
    Data = file.read()
    Data = getSequencefromfile(Data)
    return Data

def nucleotides_count(Data):
    c = 0
    g = 0
    t = 0
    a = 0
    for nucleotido in Data:
        if nucleotido =='c':
            c+=1
        elif nucleotido =='t':
            t+=1
        elif nucleotido =='g':
            g+=1
        elif nucleotido =='a':
            a+=1
    return a,c,g,t

def calculaFactorial(n):
    if n==1:
      return 1
    if n>1:
         n1 = n * calculaFactorial(n - 1)
    return n1

def create_pairbase(sequence):
    add=0
    reading_frames = []
    new_sequence = ""
    read_frame = ""
    for i in sequence:
      if i=='a':
        new_sequence += 't'
      if i=='t':
        new_sequence += 'a'
      if i=='c':
        new_sequence += 'g'
      if i=='g':
        new_sequence += 'c'
    return ''.join(reversed(new_sequence))

def create_particles(particles,sequence):
  def generate_alternate_sequence(sequence):
    def two_integernotequal(sequence):
      x = random.randint(1, len(sequence))-1
      y = random.randint(1, len(sequence))-1
      if(x==y):
        return two_integernotequal(sequence)
      else:
        return x,y
    def replace_at(cadena, idx, char):
      return cadena[:idx] + char + cadena[idx+1:]
    new_sequence = sequence
    posx,posy=two_integernotequal(sequence)
    nucleotidox = new_sequence[posx]
    nucleotidoy = new_sequence[posy]
    new_sequence = replace_at(new_sequence,posx,nucleotidoy)
    new_sequence = replace_at(new_sequence,posy,nucleotidox)
    return new_sequence,[posx,posy,nucleotidox,nucleotidoy]


  list_of_change = []
  particles_list = []
  particles_list.append(sequence)
  for particle in range(particles):
    save,lista = generate_alternate_sequence(sequence)

    if save not in particles_list:
      particles_list.append(save)
      sequence = save

  return particles_list

def make_randomsequenceordersfromsequencestring(Data,quantity):
    Generated_sequences = create_particles(quantity,Data)
    return Generated_sequences

def writeInfileGeneratedData(Data,filename):
   chains = []
   with open("GeneratedSequences\-"+filename,"r") as file:
      for line in file:
        chains.append(line.strip())
   

   print(len(chains))
   with open("GeneratedSequences\-"+filename,"a") as file:
      for i in Data:
        if(i not in chains):
            file.write(i)
            file.write("\n")

def readfileGeneratedData(filename):
  chains = []
  with open("GeneratedSequences\-"+filename,"r") as file:
    for line in file:
      chains.append(line.strip())
  return chains[len(chains)-1]

"""
geneticcode = get_nucleotidesequencefrom_gbfile("PorcineCircovirus.gb")

for i in geneticcode:
  ORFSL = getORFS_fromsequence(i,minaminoacidlength=25,include_front_ORFS=False)
  print("Virus")
  for ORFS in ORFSL:
     print(ORFS)"""
geneticcode,Organism=get_nucleotidesequence_and_metadata_from_gbfile("sequence (3).gb")
geneticcode=get_nucleotidesequencefrom_gbfile("sequence (3).gb")

for i in geneticcode:
  ORFSL = getORFS_and_metadata_fromsequence(i,minaminoacidlength=25)
  print(ORFSL)