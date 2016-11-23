#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 16:31:59 2016
python 3.5.2
BioPython 1.68

@author: ymseah
"""
from Bio import SeqIO
#from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
import re

class analyze:
    def fastq_to_fasta(self, fastq_file, fasta_file):
        """Convert fastq reads to fasta format"""
        SeqIO.convert(fastq_file, "fastq", fasta_file, "fasta")
    
    def qblast_fasta(self, fasta_file):
        """Iterate through fasta file to individually BLAST each read"""
        for read in SeqIO.parse(fasta_file, "fasta"):
            result_handle = NCBIWWW.qblast("blastn", "nt", read.format("fasta"), hitlist_size=3)
            with open("my_blast.xml", "a") as blast_results:
                blast_results.write(result_handle.read())

    def blastn_remote(self, blast_path, fasta_file_path):
        """Uses BLAST+ command line to search NCBI servers
        ./blastn -db nt -query <input_file> -out <output_file> -remote
        """
        pass
    
    def parse_blast(self, blast_out):
        """Argument: blastn output text file name
        1. Records number of queries
        2. Records top (first) blast hit to D. vulgaris, M. maripaludis, & others.
        3. Records no hits found.
        4. Output: 2 files 
            a. Number of D.v., M.m., and no hits, with list of results.
            b. Number of other hits, with list of results.
            """
        with open(blast_out) as blast_results:
            file_lines = blast_results.readlines()
        num_queries = 0
        desulfo_v_hits = []
        methano_m_hits = []
        no_hits_counter = 0
        other_hits = []
        line_counter = 0
        while line_counter < len(file_lines):
            if re.match("Query=", file_lines[line_counter]):
                num_queries += 1
            elif re.search("Sequences producing significant alignments:", file_lines[line_counter]):
                top_hit_index = line_counter + 2
                if re.search("Desulfovibrio vulgaris", file_lines[top_hit_index]):
                    desulfo_v_hits.append(file_lines[top_hit_index])
                elif re.search("Methanococcus maripaludis", file_lines[top_hit_index]):
                    methano_m_hits.append(file_lines[top_hit_index])
                else:
                    other_hits.append(file_lines[top_hit_index])                
            elif re.search("No hits found", file_lines[line_counter]):
                no_hits_counter += 1
            line_counter += 1
        print("\nDesulfovibrio vulgaris hits / Total queries: " + str(len(desulfo_v_hits)) + "/" + str(num_queries) 
              + "\nMethanococcus maripaludis hits / Total queries: " + str(len(methano_m_hits)) + "/" + str(num_queries) 
              + "\nNo hits / Total queries: " + str(no_hits_counter) + "/" + str(num_queries) 
              + "\nOther hits / Total queries: " + str(len(other_hits)) + "/" + str(num_queries) + "\n")
        output_filename_string = blast_out + "_parsed.txt"
        output_other_hits_filename_string = blast_out + "_parsed_others.txt"
        with open(output_filename_string, "w") as outfile:
            outfile.write("Desulfovibrio vulgaris hits / Total queries: " + str(len(desulfo_v_hits)) + "/" + str(num_queries) 
                  + "\nMethanococcus maripaludis hits / Total queries: " + str(len(methano_m_hits)) + "/" + str(num_queries)
                  + "\nNo hits / Total queries: " + str(no_hits_counter) + "/" + str(num_queries) + "\n\n")
            for dv in desulfo_v_hits:
                outfile.write(dv)
            for mm in methano_m_hits:
                outfile.write(mm)
        with open(output_other_hits_filename_string, "w") as outfile2:
            outfile2.write("Other hits / Total queries: " + str(len(other_hits)) + "/" + str(num_queries) + "\n\n")
            for hits in other_hits:
                outfile2.write(hits)

#fastq_to_fasta("sic_M1_pear.assembled.unmatched.fastq","sic_M1_pear.assembled.unmatched.fasta")
#parse_blast("sic_M1_pear.assembled.unmatchedBLAST.out") 
#qblast_fasta("sic_HA3.45_1.unmatched.fasta")

'''
    def check_adapter(fasta_file):
        i_transposase = {'read1': 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG', 
                         'read2': 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'}
        i_neXT_i7v1 = {'N701': 'TAAGGCGA', 'N702': 'CGTACTAG', 
                       'N703': 'AGGCAGAA', 'N704': 'TCCTGAGC', 
                       'N705': 'GGACTCCT', 'N706': 'TAGGCATG', 
                       'N707': 'CTCTCTAC', 'N708': 'CAGAGAGG',
                       'N709': 'GCTACGCT', 'N710': 'CGAGGCTG',
                       'N711': 'AAGAGGCA', 'N712': 'GTAGAGGA'}
        i_neXT_i7v2 = {'N714': 'GCTCATGA', 'N715': 'ATCTCAGG',
                       'N716': 'ACTCGCTA', 'N718': 'GGAGCTAC',
                       'N719': 'GCGTAGTA', 'N720': 'CGGAGCCT',
                       'N721': 'TACGCTGC', 'N722': 'ATGCGCAG',
                       'N723': 'TAGCGCTC', 'N724': 'ACTGAGCG',
                       'N726': 'CCTAAGAC', 'N727': 'CGATCAGT',
                       'N728': 'TGCAGCTA', 'N729': 'TCGACGTC'}
        i_neXT_i5v1_miseq = {'S501': 'TAGATCGC', 'S502': 'CTCTCTAT',
                             'S503': 'TATCCTCT', 'S504': 'AGAGTAGA',
                             'S505': 'GTAAGGAG', 'S506': 'ACTGCATA',
                             'S507': 'AAGGAGTA', 'S508': 'CTAAGCCT'}
        i_neXT_i5v2_miseq = {'S510': 'CGTCTAAT', 'S511': 'TCTCTCCG',
                             'S513': 'TCGACTAG', 'S515': 'TTCTAGCT',
                             'S516': 'CCTAGAGT', 'S517': 'GCGTAAGA',
                             'S518': 'CTATTAAG', 'S520': 'AAGGCTAT',
                             'S521': 'GAGCCTTA', 'S522': 'TTATGCGA'}
    
        for read in SeqIO.parse(fasta_file, "fasta"):
            pass
        pass
'''