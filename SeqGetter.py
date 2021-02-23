from Bio import Entrez
from Bio import SeqIO
Entrez.email = "jcolgan@luc.edu"
output = open("humanHerpesvirusCompleteSequence.txt", 'a')
handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
count =0
for feature in record.features:
        if feature.type == "CDS":
            count+=1
            feature_name = "feature.qualifiers" 
            feature_seq = feature.extract(record.seq)
            # Simple FASTA output without line wrapping:
            output.write(">" + feature_name + "\n" + str(feature_seq) + "\n")
output.close()
output=open("miniProjectLog",'a')
output.write('The HCMV genome (EF999921) has '+str(count) +' CDS.')
