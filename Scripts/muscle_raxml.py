import os
import subprocess
import sys

base_dir=sys.argv[1]+"/prod_fasta/"
os.chdir(base_dir+'seqrecords/pergene_seqrecords')

os.system("mkdir muslce_output")
for file in os.listdir("."):
     if file.endswith(".fasta"): 
         base=os.path.basename(file)
         filename=os.path.splitext(base)[0]
         file=os.path.join(base_dir+'seqrecords/pergene_seqrecords',file)
         subprocess.call(["muscle","-in",file,"-out","./muslce_output"+"/"+filename+"_"+"align"+".fasta","-diags"])