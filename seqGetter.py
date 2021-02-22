
m Bio import Entrez
from Bio import SeqIO
Entrez.email = "jcolgan@luc.edu"
output = open("Human herpesvirus 5 strain TB40/E clone TB40-BAC4, complete sequence", w)
        handle = Entrez.efetch(db="nucleotide", id="EF999921", rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        output.write(str(record.seq))
        output.close

        
