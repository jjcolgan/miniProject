from Bio import Entrez
from Bio import SeqIO
import os
import argparse
parser = argparse.ArgumentParser(description='Switch between test and real pipeline.')
parser.add_argument('test', metavar='t', type=str, nargs='+',help='a switch between test and real pipeline')
parser.add_argument('full', metavar='f', type=str, nargs='+',help='a switch between test and real pipeline')
args=parser.parse_args()



if args.full:
    #create dir for the files nd such and change to it 
    os.system('mkdir miniProject_JJ_Colgan')
    os.chdir('miniProject_JJ_Colgan')
    os.system('touch miniProjectLog')
    wget = 'wget '
    # download the data from sra using wget
    os.system(wget + 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660030/SRR5660030.1') 
    os.system(wget + 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660033/SRR5660033.1')
    os.system(wget + 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660044/SRR5660044.1')
    os.system(wget + 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR5660045/SRR5660045.1')

    fastqDump = 'fastq-dump -I --split-files '
    # Decompressing fastq files into the forward and reverse reads 
    os.system(fastqDump+ 'SRR5660030.1')
    os.system(fastqDump+ 'SRR5660033.1')
    os.system(fastqDump+ 'SRR5660044.1')
    os.system(fastqDump+ 'SRR5660045.1')
if args.test:
    os.system('mkdir miniProject_JJ_Colgan'):wq
    os.system('mv PracticeSRR5660030.1_1.fastq SRR5660030.1_1.fastq')
    os.system('mv PracticeSRR5660030.1_2.fastq SRR5660030.1_2.fastq')
    os.system('mv PracticeSRR5660033.1_1.fastq SRR5660030.1_1.fastq')
    os.system('mv PracticeSRR5660033.1_2.fastq SRR5660030.1_1.fastq')
    os.system('mv PracticeSRR5660044.1_1.fastq SRR5660030.1_1.fastq')
    os.system('mv PracticeSRR5660044.1_2.fastq SRR5660030.1_1.fastq')
    os.system('mv PracticeSRR5660045.1_1.fastq SRR5660030.1_1.fastq')
    os.system('mv PracticeSRR5660045.1_2.fastq SRR5660030.1_1.fastq')
    os.system('mv SRR5660030.1_1.fastq miniProject_JJ_Colgan')
    os.system('mv SRR5660030.1_2.fastq miniProject_JJ_Colgan')
    os.system('mv SRR5660033.1_1.fastq miniProject_JJ_Colgan')
    os.system('mv SRR5660033.1_2.fastq miniProject_JJ_Colgan')
    os.system('mv SRR5660044.1_1.fastq miniProject_JJ_Colgan')
    os.system('mv SRR5660044.1_2.fastq miniProject_JJ_Colgan')
    os.system('mv SRR5660045.1_1.fastq miniProject_JJ_Colgan')
    os.system('mv SRR5660045.1_2.fastq miniProject_JJ_Colgan')
    os.chdir('miniProject_JJ_Colgan')
    os.system('touch miniProjectLog')

# retriving the transciptome index for HCMV
# count is used to count the CDS
# basically did an entrez fetch for the id EFF999921, build a record, and read using seqIO. Then used a for loop which looked for CDS in features. IF there was a CDS, it was turned into a fasta
# record using > feature name, and count was incremented. This record is added to the output called humanHerpesvirusCompleteSequence.txt. After the loop is ended, the log is opened and the
# count recorded. 
Entrez.email = "jcolgan@luc.edu"
output = open("humanHerpesvirusCompleteSequence.txt", 'a')
handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
count =0
for feature in record.features:
    if feature.type == "CDS":
        count+=1
        feature_name = record.description
        feature_seq = feature.extract(record.seq)
        # Simple FASTA output without line wrapping:
        output.write(">" + feature_name + "\n" + str(feature_seq) + "\n")
output.close()
output=open('miniProjectLog','a')
output.write('The HCMV genome (EF999921) has '+str(count) +' CDS.\n')
output.close()

#creating index of transcriptome from humanHerpesvirusCompleteSequence

os.system('time kallisto index --make-unique -i humanHerpesvirusIndex.idx humanHerpesvirusCompleteSequence.txt')

#Quantification of each sample
kallisto = "time kallisto quant -i humanHerpesvirusIndex.idx -b 3 -t 4 -o "
fastq1 = '_1.fastq '
fastq2 = '_2.fastq'
os.system(kallisto + 'SRR5660030.1quant ' +'SRR5660030.1'+ fastq1 + 'SRR5660030.1' + fastq2)
os.system(kallisto + 'SRR5660033.1quant ' +'SRR5660033.1' + fastq1 + 'SRR5660033.1'+ fastq2)
os.system(kallisto + 'SRR5660044.1quant '+ 'SRR5660044.1'+ fastq1 + 'SRR5660044.1'+ fastq2)
os.system(kallisto + 'SRR5660045.1quant '+ 'SRR5660045.1'+ fastq1 + 'SRR5660045.1'+ fastq2)
#Write results to log / format log
output = output=open('miniProjectLog','a')
output.close()
# building table for sleuth
os.system('touch table.txt')
output = open('table.txt','w')
output.write('sample\t condition\t path\n')
output.write('SRR5660030.1\t 2dpi\t SRR5660030.1quant\n')
output.write('SRR5660033.1\t 6dpi\t SRR5660033.1quant\n')
output.write('SRR5660044.1\t 2dpi\t SRR5660044.1quant\n')
output.write('SRR5660045.1\t 6dpi\t SRR5660045.1quant\n')
output.close()
#buiding r code to run sleuth
os.system('touch rScript.R')
output=open('rScript.R','w')
output.write('library(sleuth)\n')
output.write('library(dplyr)\n')
output.write('stab <- read.table("table.txt",header=TRUE,stringsAsFactors=FALSE)\n')
output.write('so <- sleuth_prep(stab)\n')
output.write("so <- sleuth_fit(so, ~condition, 'full')\n")
output.write("so <- sleuth_fit(so, ~1, 'reduced')\n")
output.write("so <- sleuth_lrt(so, 'reduced', 'full')\n")
output.write("sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)\n")
output.write("sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)\n")
output.write("write.table(sleuth_significant, file='rResults',quote = FALSE,row.names = FALSE,sep ='\t')\n")
output.close()
#running sleuth 
os.system('Rscript rScript.R')
#reading results of r into the log read all lines into a list, then split each element by the tab regex. We know the first 4 items of the table are what we want, so write the first 4
# from the split line into the file 
input = open('rResults','r')
lines=[]
lines=input.readlines()
output=open('miniProjectLog','a')
temp=[]
for line in lines:
   temp = line.split('\t')
   output.write(temp[0]+'\t' +temp[1]+'\t' + temp[2] +'\t' +temp[3]+'\n')
output.close()
input.close()


#buidling bowtie2 index
os.system('bowtie2-build humanHerpesvirusCompleteSequence.txt HCMV')

##filtering reads that do not map to ind)ex
os.system('bowtie2 -x HCMV -1 SRR5660030.1'+fastq1 + '-2 SRR5660030.1' + fastq2 + ' -S SRR5660030.1.sam --al-conc SRR5660030.1_mapped_.fq')
os.system('bowtie2 -x HCMV -1 SRR5660033.1'+fastq1 + '-2 SRR5660033.1' + fastq2 + ' -S SRR5660033.1.sam --al-conc SRR5660033.1_mapped_.fq')
os.system('bowtie2 -x HCMV -1 SRR5660044.1'+fastq1 + '-2 SRR5660044.1' + fastq2 + ' -S SRR5660044.1.sam --al-conc SRR5660044.1_mapped_.fq')
os.system('bowtie2 -x HCMV -1 SRR5660045.1'+fastq1 + '-2 SRR5660045.1' + fastq2 + ' -S SRR5660045.1.sam --al-conc SRR5660045.1_mapped_.fq')
a = 0
b = 0
#counting the number of reads before and after, count is determined by the number of lines starting with @ as those are the beginning of fastq seqs. Both the forward and reverse
#should have the same number mapped, so only count one
output=open('miniProjectLog','a')
input = open('SRR5660030.1_1.fastq','r')
lines = input.readlines()

for line in lines:
    if line.startswith('@'):
        a +=1
lines.clear()
input.close()
input = open('SRR5660030.1_mapped_.1.fq','r')
lines = input.readlines()
for line in lines:
    if line.startswith('@'):
        b+=1
input.close()
lines.clear()
output.write('Donor 1 (2dpi) had ' +str(a)+ ' read pairs before Bowtie2 filtering and '+ str(b)+' read pairs after.\n')
a=0
b=0
input = open('SRR5660033.1_1.fastq','r')
lines = input.readlines()
for line in lines:
    if line.startswith('@'):
        a +=+1
input.close()
lines.clear()
input = open('SRR5660033.1_mapped_.1.fq','r')
lines = input.readlines()
for line in lines:
    if line.startswith('@'):
        b+=1
input.close()
lines.clear()
output.write('Donor 1 (6dpi) had ' +str(a)+ ' read pairs before Bowtie2 filtering and '+ str(b)+' read pairs after.\n')
input = open('SRR5660044.1_1.fastq','r')
lines = input.readlines()
a=0
b=0
for line in lines:
    if line.startswith('@'):
        a +=1
input.close()
lines.clear()
input = open('SRR5660044.1_mapped_.1.fq','r')
lines = input.readlines()
for line in lines:
    if line.startswith('@'):
        b+=1
input.close()
lines.clear()
output.write('Donor 2 (2dpi) had ' +str(a)+ ' read pairs before Bowtie2 filtering and '+str(b)+' read pairs after.\n')
a=0
b=0
input = open('SRR5660044.1_1.fastq','r')
lines = input.readlines()
for line in lines:
    if line.startswith('@'):
        a +=1
input.close()
lines.clear()
input = open('SRR5660045.1_mapped_.1.fq','r')
lines = input.readlines()
for line in lines:
    if line.startswith('@'):
        b+=1
input.close()
output.write('Donor 2 (2dpi) had ' +str(a)+ ' read pairs before Bowtie2 filtering and '+ str(b)+' read pairs after.\n')
output.close()
lines.clear()
#assembling with spades
os.system('spades -k 55,77,99,127 -t 2 --only-assembler -1 SRR5660030.1_mapped_.1.fq -2 SRR5660030.1_mapped_.2.fq -1 SRR5660033.1_mapped_.1.fq -2 '
          'SRR5660033.1_mapped_.2.fq -1 SRR5660045.1_mapped_.1.fq -2 SRR5660045.1_mapped_.2.fq -1 SRR5660044.1_mapped_.1.fq -2 SRR5660044.1_mapped_.2.fq -o assembly/')
#writing spades command to log 
output=open('miniProjectLog','a')
output.write('spades -k 55,77,99,127 -t 2 --only-assembler -1 SRR5660030.1_mapped_.1.fq -2 SRR5660030.1_mapped_.2.fq -1 SRR5660033.1_mapped_.1.fq -2 '
          'SRR5660033.1_mapped_.2.fq -1 SRR5660045.1_mapped_.1.fq -2 SRR5660045.1_mapped_.2.fq -1 SRR5660044.1_mapped_.1.fq -2 SRR5660044.1_mapped_.2.fq -o assembly/\n')
output.close()
#counting fasta seqs that are longer than 1000 bp using a to count the number of seqs over 1000, if a seq is over 1000, its length is added to b which is used to keep track of the
#assembly length. Max len is used to keep track of the longest length seq there is in the assembly. If a seq is encounter whos length is longer, its length is set to maxLen and its
# seq made maxlenstr. this is then wrote to the file longestContig.fasta, and forced into fasta format with the title, >longboi\n followed by the seq. This is then used to perform a blast
# from which a tab demlimited file is made of the output by splitting along the commas. From this the first 10 are taken and written to log using a for loop/ counter. 
input = open('assembly/contigs.fasta','r')
a=0
b=0
maxLen =0
maxlenstr=''
for record in SeqIO.parse(input, 'fasta'):
    if len(record.seq) >  1000:
        a +=1
        b =len(record.seq)+b
        if len(record.seq) >  maxLen:
            maxLen = len(record.seq)
            maxlenstr = str(record.seq)
input.close()
output=open('miniProjectLog','a')
output.write('There are ' + str(a)+ " contigs > 1000 bp in the assemby.\n")
output.write('There are '+ str(b) + ' bp in the assembly.')
os.system('touch longestContig.fasta')
output= open('longestContig.fasta','w')
output.write('>longboi\n')
output.write(maxlenstr)
output.close()
os.system('blastn -query longestContig.fasta -db viraldb -out blastnResults -outfmt='
          '"10 sacc pident length qstart qend sstart send bitscore evalue stitle"')
input = open('blastnResults', 'r')
output = open('miniProjectLog','a')
output.write('sacc\t pident\t length\t qstart\t qend\t sstart\t send\t bitscore\t evalue\t stitle\n')
lines = input.readlines()
counter = 0
for line in lines:
    temp = line.split(',')
    output.write(temp[0]+'\t'+temp[1]+'\t'+temp[2]+'\t'+temp[4]+'\t'+temp[5]+'\t'+temp[6]+temp[7]+'\t'+temp[8]+temp[9])
    counter+=1
    if counter == 9:
        output.close()
        break 

    




