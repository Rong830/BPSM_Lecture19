#!/usr/local/bin/python3

import Bio
from Bio import Entrez, SeqIO
import numpy as np
import pandas as pd
import os, glob
import subprocess
import datetime

# Settings
Entrez.email = "s2160628@ed.ac.uk"
Entrez.api_key=subprocess.check_output("echo $NCBI_API_KEY", shell=True).rstrip().decode()

'''
Questions
1. How many complete COX1 protein records are there for mammals?
2. What is their average length (the proteins that is, not the mammals!)?
3. Write a function that will answer the question for any gene name and any taxonomic group.
4. What else could you sensibly include in the function (in terms of retrieving possibly useful data)?
'''


## 1. How many complete COX1 protein records are there for mammals?
protein_result_mammals = Entrez.read(Entrez.esearch(db="protein", term="Mammalia COX1 complete", retmax="20"))
print('The number of complete COX1 protein records is: ', protein_result_mammals['Count'])


## 2. What is their average length (the proteins that is, not the mammals!)?
loopcounter = total_length = 0

for accession in protein_result_mammals['IdList']:
    loopcounter += 1
    gb_file = Entrez.efetch(db="protein",id=accession,rettype="gb")
    record = SeqIO.read(gb_file, "genbank")
    total_length =  total_length + len(record.seq)
    mean_length = int(total_length/loopcounter)
    print(record.id+"\t"+record.description+"\t"+str(len(record.seq))+"\t"+record.seq[0:10]+"...")

print("The mean length was "+str(mean_length)+" amino acids.")


## 3. Write a function that will answer the question for any gene name and any taxonomic group.
search_results = []
def get_average_length(taxonomy, gene, howmany=10):
    from Bio import Entrez, SeqIO
    import os, subprocess
    Entrez.email = "s2160628@ed.ac.uk"
    Entrez.api_key=subprocess.check_output("echo $NCBI_API_KEY", shell=True).rstrip().decode()
    search_term = taxonomy + " " + gene + " complete"
# Set up output file, named after search used (spaces removed!)
    search_output = open(search_term.replace(" ","_")+"_outputs.txt","w")
    mysearch = Entrez.esearch(db="protein", term=search_term, retmax=howmany)
    result = Entrez.read(mysearch)
    loopcounter = total_length = 0
# Extract info from the results of the search
    for accession in result['IdList']:
        loopcounter += 1
        gb_file = Entrez.efetch(db="protein",id=accession,rettype="gb")
        record = SeqIO.read(gb_file, "genbank")
        total_length =  total_length + len(record.seq)
# Add the results of the search to search_results
        search_results.append([search_term,record.id,record.description,len(record.seq),record.seq])
# Show the results (trimmed sequence) of the search to screen as we go
        print(record.id+"\t"+record.description+"\t"+str(len(record.seq))+"\t"+record.seq[0:50]+"...")
# Write results to search output file
        search_output.write(record.id+"\t"+record.description+"\t"+str(len(record.seq))+"\t"+str(record.seq)+"\n")
# Interim print statment for mean length
    mean_length = int(total_length/loopcounter)
    close(search_output)    
    return print(("\nThe mean length was "+str(mean_length)+" amino acids.\n"))

def get_average_length2(taxonomy, gene, startat=0,howmany=100):
    from Bio import Entrez, SeqIO
    Entrez.email = "s123456@ed.ac.uk"
    Entrez.api_key=subprocess.check_output("echo $NCBI_API_KEY", shell=True).rstrip().decode()
    search_term = taxonomy + " " + gene + " complete"
    search_output = open(search_term.replace(" ","_")+"_outputs.txt","w")
    mysearch = Entrez.esearch(db="protein", term=search_term, retstart = startat, retmax=howmany)
    result = Entrez.read(mysearch)
    print("Search done,", result['Count'],"found, starting retrieval of",howmany,"starting at",startat)
    loopcounter = poorseqs = total_length = 0
    for accession in result['IdList']:
      loopcounter += 1
      genbank = Entrez.efetch(db="protein",id=accession,rettype="gb")
      print('\r' + 'Retrieving sequence ' + str(loopcounter+startat), end="")
      record = SeqIO.read(genbank, "genbank")
      Xaa = str(record.seq).count("X")
      if Xaa > 5 :
         print(" Seq",loopcounter+startat,"contains",Xaa,"unknown amino acids:",int(100*Xaa/len(record.seq)), "percent of total!")
         poorseqs += 1
      else :
         total_length =  total_length + len(record.seq)
         search_results.append([search_term,record.id,record.description,len(record.seq),record.seq])
         # print(record.id+"\t"+record.description+"\t"+str(len(record.seq))+"\t"+record.seq[0:50]+"...")
         search_output.write(record.id+"\t"+record.description+"\t"+str(len(record.seq))+"\t"+str(record.seq)+"\n")
    goodseqs=loopcounter-poorseqs
    mean_length = int(total_length/goodseqs)
    close(search_output)
    return print(("\nThe mean length of the "+str(goodseqs)+" high quality seqs was "+str(mean_length)+" amino acids.\n"))

# for loop for efetching
for from_here in batch_starts :
    print(datetime.datetime.now().strftime('## %d %b %Y, %H:%M:%S'))
    to_here = from_here + batch_size
    print("Getting the sequences",from_here+1,"to",to_here)
    all_fastas.write(Entrez.efetch(db="protein", \
    id=result['IdList'][from_here:to_here], \
    retmode="text",rettype="fasta").read())

print("\n",datetime.datetime.now().strftime('## %d %b %Y, %H:%M:%S'))