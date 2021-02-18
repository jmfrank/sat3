# -*- coding: utf-8 -*-
"""
Spyder Editor

"""

from Bio import SeqIO
import os

input_file = "example.gb" #Your GenBank file locataion. e.g C:\\Sequences\\my_genbank.gb
output_file_name = "Output.fasta" #The name out your fasta output
accession_numbers = [line.strip() for line in open('example.gb')] #the same as your input file, defines the headers for each sequence

if not os.path.exists(output_file_name): #checks for a pre-existing file with the same name as the output
    for rec in SeqIO.parse(input_file, "gb"): #calls the record for the genbank file and SeqIO (BioPython module) to parse it
        acc = rec.annotations['accessions'][0] #Defines your accession numbers
        organism = rec.annotations['organism'] #defines your organism ID
        tax_line = ("| ").join(rec.annotations['taxonomy']) #defines your taxonomy and seperates entries with a |, remove this line, the 'tax_line', and the {2} in your save for a simpler output
        for feature in rec.features: #looks for features in the genbank
            for key, val in feature.qualifiers.items(): #looks for val in the feature qualifiers
                if any("protein" in s for s in val): #Finds all the CDS which contain the word "protein" in the qualifiers. Change to 'if "Name" in val:' for protein called "name" exactly
                    with open(output_file_name, "a") as ofile: #opens the output file and "a" designates it for appending
                        ofile.write(">{0}| {1}| {2}| {3}| {4}| \n{5}\n\n".format(acc, organism, tax_line, feature.qualifiers['protein_id'][0], val[0], feature.qualifiers['translation'][0])) #Writes my FASTA format sequences to the output file
                        
else:
    print ("The output file already seem to exist in the current working directory {0}. Please change the name of the output file".format(os.getcwd())) #error code, so you don't overwirite your files