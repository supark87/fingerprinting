import os
import subprocess
import sys

script=sys.argv[0]
base_dir=sys.argv[1]
coregenome_dir=sys.argv[1]+"/../coregene/"

os.chdir(base_dir)
def MakeBlastDB():  
     os.system("mkdir output")
     for genome in os.listdir('.'):
         genome=os.path.join(base_dir,genome)
         base=os.path.basename(genome)
         basename=base.split(".")[0]
         os.system("module load ncbi-blast+/LATEST") 
         subprocess.call(["makeblastdb","-in",genome,"-dbtype","nucl","-out","./output"+"/"+basename])
#MakeBlastDB()

child_processes=[]
def RunBlast(): 
     os.system("mkdir blastoutput")
     for query in os.listdir(coregenome_dir):
         if query.endswith(".fasta"):
             query=os.path.join(coregenome_dir,query)
             baseq=os.path.basename(query)
             filename =os.path.splitext(baseq)[0] 
             for database in os.listdir(base_dir+"/output"):
                 database=os.path.join(base_dir+"/output",database)
                 basedb=os.path.basename(database)
                 print(basedb)
                 dbname=basedb.split(".")[0]
                 databasename =os.path.join(base_dir+"/output",basedb.split(".")[0])
                 p=subprocess.Popen(["blastn","-query",query,"-db",databasename,"-outfmt","6 qseqid sseqid pident qlen qstart qend sstart send","-out","./blastoutput"+"/"+filename+"_"+dbname+".blast"])
                 child_processes.append(p)
                 for cp in child_processes:
                     cp.wait()
#RunBlast()
#print("blast is done")


os.chdir(base_dir+"/blastoutput")

def filter():
     os.system("mkdir sorted_blast_pair") 
     for blastresult in os.listdir('.'):
         if blastresult.endswith(".blast"):
             genomename=os.path.basename(blastresult)
             blastresult=open(blastresult)
             for line in blastresult:
                 try:
                     gene={}
                     line = line.split( )
                     qseqid=line[0]
                     sseqid=line[1]
                     pident=float(line[2])
                     qlength=float(line[3])
                     qstart=float(line[4])
                     qend=float(line[5])
                     sstart=float(line[6])
                     sstart=float(line[6])
                     send=float(line[7])
                     if (pident>98) & (((qend-qstart+1)/qlength)>0.99) :
                         gene[qseqid]=sseqid
                         sseqname=sseqid.split("_")[0]+"_"+sseqid.split("_")[1]+"_"+sseqid.split("_")[2]
                         for key in gene:
                             name=str(key).split("_")[0]+"_"+str(key).split("_")[1]
                             #print(name)
                             with open("./sorted_blast_pair"+"/"+name+"_"+sseqname+".pair","w") as ofile:
                                  ofile.write(key+"\t"+gene.get(key))
                                  ofile.close 
                 except IOError:
                     print("no input")
             blastresult.close() 
filter()
#print("Filtering blast result is done")


####GetSequence#####

os.chdir(base_dir)
os.system("mkdir seqrecords")
def Parse(filename,seqs):
    file = open(filename) 
    seqs={}
    name = ''
    for line in file:
        line = line.rstrip()  
        if line.startswith('>'):
            name=line.replace('>',"")
            seqs[name] = ''
        else:
            seqs[name] = seqs[name] + line
            
    file.close
    return seqs

seqs={}    
for genome in os.listdir('.'):

     if genome.endswith(".fasta"):
         seqs=dict(seqs,**Parse(genome,seqs))
            
for file in os.listdir(base_dir+'/blastoutput/sorted_blast_pair'):
     genomename=file.split("_")[2]+"_"+file.split("_")[3]+"_"+file.split("_")[4]
     coregenename=file.split("_")[0]+"_"+file.split("_")[1]
     file=open(os.path.join(base_dir+'/blastoutput/sorted_blast_pair',file))
     for line in file:
         line=line.rstrip()
         genename=line.split("\t")[1]
         #coregenename=line.split("\t")[0]
         for key in seqs:   
             if key.find(str(genename))!= -1: 
                 #print(genename)
                 with open("./seqrecords"+"/"+coregenename+"_"+genename+"_"+genomename+".fasta","w") as ofile:
                      ofile.write(">"+coregenename+"_"+genename+"_"+genomename+"\n"+seqs.get(key))
                      ofile.close()
     file.close()
    
print("Getting sequences are done")     		

os.chdir(base_dir+'/seqrecords')
os.system('mkdir pergene_seqrecords')

for (a,b,files) in os.walk("."):
    for seqrecord in files:
        if seqrecord.endswith("fasta"):
            gene=seqrecord.split("_")[0]+"_"+seqrecord.split("_")[1]
            seq=open(os.path.join(base_dir+"/seqrecords",seqrecord),"r")
        #seq=open(seqrecord,"r")
            for seqline in seq:
                seqline=seqline.rstrip()
                with open(gene+"_"+"unaligned"+".fasta","a") as pfile:
                    pfile.write(seqline+"\n")
                    pfile.close 
            seq.close()
os.system("mv *unaligned.fasta ./pergene_seqrecords")

print("Sequences are sorted by each locus")


#####PRESENCE/ABSENCE#####
os.chdir(base_dir)
filelist1=os.listdir(base_dir+"/blastoutput")
filelist2=os.listdir(base_dir+"/blastoutput/sorted_blast_pair")
sys.stdout=open('test','a')
for before in filelist1:
     before_name=before.split("_")[0]+"_"+before.split("_")[1]+"_"+before.split("_")[2]
     if before.endswith(".blast"):
         base=before.split("blast")[0]
         coregenename=base.split("_")[0]+"_"+base.split("_")[1]
         #print(coregenename)
         genomename=base.split("_")[2]
         #print(genomename)
         for after in filelist2:
             after1=after.split("_")[0]+"_"+after.split("_")[1]+"_"+after.split("_")[2]
             #print(after1)
             after2=after.split("_")[2]+"_"+after.split("_")[3]+"_"+after.split("_")[4]
			 #print(after2)
             if str(after).find(str(coregenename))!= -1:
                 sys.stdout.write(coregenename+"\t"+after2+"\t"+"yes"+"\n")      
             else:
                 sys.stdout.write(coregenename+"\t"+after2+"\t"+"no"+"\n")
         
sys.stdout.close()

base_dir=sys.argv[1]
os.chdir(base_dir+'/'+'/seqrecords/pergene_seqrecords')
os.system("mkdir muslce_output")
for file in os.listdir("."):
     if file.endswith(".fasta"): 
         base=os.path.basename(file)
         filename=os.path.splitext(base)[0]
         file=os.path.join(base_dir+'seqrecords/pergene_seqrecords',file)
         subprocess.call(["muscle","-in",file,"-out","./muslce_output"+"/"+filename+"_"+"align"+".fasta","-diags"])
