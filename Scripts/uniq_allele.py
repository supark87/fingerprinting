from Bio import AlignIO
import os
import sys

base_dir=sys.argv[1]
if not os.path.exists('allele_subtype'):
    os.makedirs('allele_subtype')
for file in os.listdir(base_dir):
     if file.endswith(".fasta"):
         corename=file.split("_")[0]
         count=1
         alignment=AlignIO.read(file,"fasta")
         dic1={}
         for record in alignment :
              dic1[str(record.seq)]=1
         for key1 in dic1:
              with open(file+"_"+"uniq.fasta",'a') as ofile:
                  ofile.write(">"+corename+"_"+str(count)+"\n"+key1+"\n")
              count=count+1
command="mv *_uniq.fasta ./allele_subtype/"
os.system(command)