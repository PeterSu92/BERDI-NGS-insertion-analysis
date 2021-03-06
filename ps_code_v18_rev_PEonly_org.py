#Import all libraries of interest
import gzip 
import logging
import os
import re
import sys
import time
import uuid
import random
import numpy as np
import pandas as pd
import csv
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Emboss.Applications import NeedleCommandline
import sys, getopt, os
#Code subset up through filter_pe_mismatch only
from util.regex_compile.functions import compile_res, filter_seqs
# from util.insertion_site_calculation import insertion_chunks, insertion_site_freq
from util.pe_matching.pe_copying import get_coords, get_sense, gen_copied_seq_function, get_copied_seq, filter_pe_mismatch, score_cutoff_by_length,quality_filter_single 
from util.filtering.functions import quality_filter#,cull_alignments,alignment_filter 

def load_ngs_file(fpath,ftype='fastq'):
    """
    Load a .fastq file to a SeqIO iterator, un gzip if necessary, as files are
    often returned to us as .gzip files.
    
    inputs:
        fpath - str, path to the file (either gzip or fastq)
        ftype = str, in this case always 'fastq' since NGS files are all FASTQ
        
    outputs:
        f_iter - generator object that doesn't actually load the data until you
        call upon it - saves time until you actually need the data
    """
    if fpath.endswith('.gz'):
        seq_f = gzip.open(fpath,'rb')
        # The above line opens a gzip-compressed file in binary or text mode
        # the 'rt' indicates it is opened in text mode; 'rb' indicates binary mode
    elif fpath.endswith('.fastq'):
        seq_f = open(fpath,'rb')
        # If it finds any fastq files, it simply opens them in text format (or binary)
    else:
        raise ValueError('File does not end in .gz or .fastq; confirm file type')
    f_iter = SeqIO.parse(seq_f,ftype)
    return f_iter
    
    #We also need a regex function (efficient ways of searching for string) for compiling 
    #This is essentially faster since we will compile the sequences first
    

            
#Next we want to filter through all of the fastq file to only get the sequences
#that have the corresponding required sequences in f_res
            

            
#Next we will start filtering based on the sequences desired in compile_res

#This is the part that gets run
def filter_sample(f_name,pe_name,template,f_filt_seqs,r_filt_seqs):

        f_res = compile_res(f_filt_seqs)
        pe_res = compile_res(r_filt_seqs)

        # Load FASTQ files as generator
        f_seqs = load_ngs_file(f_name)
        pe_seqs = load_ngs_file(pe_name)
        # Create lists of SeqRecord objects from the generators
        #Sometimes, the NGS reads have ambiguous bases or improperly sequenced bases, simply labled 'N'. We cannot use these reads and so they are removed here
        f_seqs1 = [s for s in f_seqs if 'N' not in str(s.seq)]
        pe_seqs1 = [s for s in pe_seqs if 'N' not in str(s.seq)]
        print(str(len(f_seqs1))+' forward reads and '+str(len(pe_seqs1))+' paired end reads initially')  
        #Filter for quality
        # These sequences with "2" at the end will be filtered for the sequences
        # corresponding to the MBP primer, ZFP primer, and the transposon scar
        f_seqs2 = []
        pe_seqs2 = []
        for regex in f_res:
            f_seqs2.append(filter_seqs(f_seqs1,regex))
                
        # As for the forward reads, repeat for the paired-end reads
        for regex in pe_res:
            pe_seqs2.append(filter_seqs(pe_seqs1,regex))
        print('Forward reads:'+str(len(f_seqs2[0]))+',',str(len(f_seqs2[1]))+',',str(len(f_seqs2[2])))
        print('PE reads:'+str(len(pe_seqs2[0]))+',',str(len(pe_seqs2[1]))+',',str(len(pe_seqs2[2])))
        #At this point, both f_seqs2 and pe_seqs2 will be a list of lists.

        #Depending on the direction the code searches in, a particular sequence is required for the both the forward and paired end reads. 
        f_seqs3 = f_seqs2[0]
        # Repeat for paired-end reads
        pe_seqs3 = pe_seqs2[2]
        # f_seqs3 = f_seqs3[0:1000]
        # pe_seqs3=pe_seqs3[0:1000]
       	f_seqs =[]
       	pe_seqs = []
        print(str(len(f_seqs3))+' forward reads have the sequence of interest (MBP forward primer)')
        print(str(len(pe_seqs3))+' paired-end reads have the sequence of interest (transposon scar)')
        # Now that only sequences containing the CS and the TR have been filtered for, respectively
        # the paired-end matching can occur
        f_seqs1 = []
        pe_seqs1 = []
        f_seqs2 = []
        pe_seqs2 = []
        s1 = filter_pe_mismatch(f_seqs3,pe_seqs3,gen_copied_seq_function(f_res),f_filt_seqs[2]) #right now using the scar
        seqs = s1[0]
        with open('matched_seq_PE.fa','w') as sh: #create temporary seq file, hopefully re-written each time
                #temp_f_seq = copied
             SeqIO.write(seqs,sh,'fastq')
        read_len_postalign_list = s1[1]
        print(str(len(seqs))+' forward reads have a paired-end match')

        with open('read_lengths_PE.csv','w') as file:
            writer = csv.writer(file)
            writer.writerow(["pe_append","phred_len"])
            writer.writerows(s1[1])
            file.close()

        seqs = quality_filter(seqs,q_cutoff=20)
        print(str(len(seqs))+' forward reads survived the Phred score quality filter')

        return seqs



def main(argv):
    template_file = 'temptemplate.fa'
    template = str(list(SeqIO.parse(template_file,'fasta'))[0].seq)
    f_filt_seqs_file1 = ''
    f_name = ''
    pe_name = ''
    outp_name = ''

    try:
        opts, args = getopt.getopt(argv,"hc:f:p:o:",["csvfile=","ffile=","pefile=", "outpfile="])
    except getopt.GetoptError:
        print ('process_fileset.py -c <csvfile> -f <f_file> -p <pe_file> -o <out_file>')
        sys.exit(2)
    print(opts)
    for opt, arg in opts:
        if opt == '-h':
            print ('process_fileset.py -c <csvfile> -f <f_file> -p <pe_file> -o <out_file>')
            sys.exit()
        elif opt in ("-c", "--csvfile"):
            f_filt_seqs_file1 = arg
        elif opt in ("-f", "--ffile"):
            f_name = arg
        elif opt in ("-p", "--pefile"):
            pe_name = arg
        elif opt in ("-o", "--outpfile"):
            outp_name = arg

    print ('CSV File is"', f_filt_seqs_file1)
    print ('f_name file is "', f_name)
    print ('pe_name file is "', pe_name)
    print ('CSV output file is', outp_name)

    output_file_prefix = os.path.basename(os.path.splitext(f_filt_seqs_file1)[0])
    # this line strips something like './filter_sequences/rxn1_828_829_F.csv' to a simple string 'rxn1_828_829_F'
    # and that get used as a PREFIX of a file name for the figure file and result csv file below.

    f_df = pd.read_csv(f_filt_seqs_file1)
    f_filt_seqs = f_df['sequence'].tolist()
    print ('f_filt_seqs ', f_filt_seqs)
    # f_filt_seqs  ['CAACGATTCATACATAGCTAAAAGGTACC', 'CTCTAGGGCTAGCTCTAGCCAT', 'TGCGGCCGCA']
    r_filt_seqs = []
    #Generate reverse compliment for r_filt_seqs
    r_filt_seqs = [str(Seq(seq).reverse_complement()) for seq in f_filt_seqs]
    print ('r_filt_seqs ', r_filt_seqs)
    final_sequences = filter_sample(f_name,pe_name,template,f_filt_seqs,r_filt_seqs)
    print(str(len(final_sequences))+' forward reads survived the final filtering')

if __name__ == "__main__":
   main(sys.argv[1:])