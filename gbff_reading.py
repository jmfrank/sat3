
import Bio
from Bio import SeqIO

# The provided gbff file does not appear to have gene annotations... useless!
gbff_file = "/media/ngs/data/genomes/HG002_CCS_canu_paternal/GCA_004796285.1_HG002_CCS_canu_paternal_1.0_genomic.gbff"

records = SeqIO.parse(gbff_file, 'genbank')
output_file_name = "Output.fasta" #The name out your fasta output
for rec in records:
    acc = rec.annotations['accessions'][0]  # Defines your accession numbers
    organism = rec.annotations['organism']  # defines your organism ID
    tax_line = ("| ").join(rec.annotations['taxonomy'])  # defines your taxonomy and seperates entries with a |, remove this line, the 'tax_line', and the {2} in your save for a simpler output
    for feature in rec.features:  # looks for features in the genbank
        print(feature)
        for key, val in feature.qualifiers.items():  # looks for val in the feature qualifiers
            if any("gene" in s for s in val):  # Finds all the CDS which contain the word "protein" in the qualifiers. Change to 'if "Name" in val:' for protein called "name" exactly
                print(rec)
                input(' ')
                with open(output_file_name, "a") as ofile:  # opens the output file and "a" designates it for appending
                    ofile.write(">{0}| {1}| {2}| {3}| {4}| \n{5}\n\n".format(acc, organism, tax_line, feature.qualifiers['gene'][0], val[0], feature.qualifiers['seq'][0]))  # Writes my FASTA format sequences to the output file