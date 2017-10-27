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
from util.insertion_site_calculation.insert_fn import insertion_chunks, insertion_site_freq,cull_alignments,alignment_filter 
from util.pe_matching.pe_copying import get_coords, get_sense, gen_copied_seq_function, get_copied_seq, filter_pe_mismatch, score_cutoff_by_length,quality_filter_single 
from util.filtering.functions import quality_filter 

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
        #Reset these lists to try and save memory/computing time
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

        # with open('read_lengths_PE.csv','w') as file:
        #     writer = csv.writer(file)
        #     writer.writerow(["pe_append","phred_len"])
        #     writer.writerows(s1[1])
        #     file.close()

        seqs = quality_filter(seqs,q_cutoff=20)
        print(str(len(seqs))+' forward reads survived the Phred score quality filter')
        bin1 = [] #bin1 is all sequences under 50 bases
        bin2 = [] #bin2 is reads between 50 and 100 bases
        bin3 = [] #bin3 is reads between 100 and 150 bases
        bin4 = [] #bin4 is reads between 150 and 200 bases
        bin5 = [] #bin5 is reads between 200 and 250 bases
        bin6 = [] #bin6 is reads between 250 and 300 bases
        bin7 = [] #bin7 contains between 300 and 350 bases
        # Create a list of all of the bins for iterative purposes
        big_bin = [bin1,bin2,bin3,bin4,bin5,bin6,bin7]
        # Cutoff scores for each bin
        bin_scores = [[42,251],[205,501],[446,751],[687,1001],[928,1251],[1100,1500],[1400,1750]]
        # Put sequences into bins based on their length        
        for s in seqs:
                
                if (len(str(s.seq).strip('-')) >= 12) and (len(str(s.seq).strip('-')) < 50): #with the primer sequences included, nothing should be below 12 bp
                    bin1.append(s)
                elif ((len(str(s.seq).strip('-')) >= 50) and (len(str(s.seq).strip('-')) < 100)):        
                    bin2.append(s)
                elif ((len(str(s.seq).strip('-')) >= 100) and (len(str(s.seq).strip('-')) < 150)):
                    bin3.append(s)                    
                elif ((len(str(s.seq).strip('-')) >= 150) and (len(str(s.seq).strip('-')) < 200)):
                    bin4.append(s)            
                elif ((len(str(s.seq).strip('-')) >= 200) and (len(str(s.seq).strip('-')) < 250)):
                    bin5.append(s)
                elif ((len(str(s.seq).strip('-')) >= 250) and (len(str(s.seq).strip('-')) < 300)):
                    bin6.append(s)
                elif ((len(str(s.seq).strip('-')) >= 300) and (len(str(s.seq).strip('-')) < 350)):
                    bin7.append(s)
        newSeqs = []
        # Run alignment with score cutoffs based on read length    
        for i in range(len(big_bin)):
                if len(big_bin[i]) == 0:
                    print ('Bin'+str(i+1)+' has no reads')
                else:
                    f = []
                    lo_cutoff = bin_scores[i][0]
                    hi_cutoff = bin_scores[i][1]
                    print('Sequences in bin '+str(i+1)+' before alignment: '+str(len(big_bin[i])))
                    f = alignment_filter(big_bin[i],template, lo_cutoff, hi_cutoff, "_bin"+str(i+1), gapopen = 10, gapextend = 0.5)        
                    print('Sequences in bin '+str(i+1)+' after alignment'+':', len(f))
                    newSeqs.append(f)
        final_seqs = []
        #this next command just takes each item from each of the sublists in newSeqs and dumps them into one new list
        final_seqs = [item for sublist in newSeqs for item in sublist]
        
        return final_seqs



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
    window1 = [1,5,9,13]
    window2 = [2,6,10,14]
    window3 = [3,7,11,15]
    window4 = [4,8,12,16]
    base_add = [-29,304,308,816] #since the template is broken into three windows, need to add the number of bases to arrive at the correct insertion
    num_regex = re.search('rxn',str(f_filt_seqs_file1))
    num_regex.span()[1]
    rxn_num = int(f_filt_seqs_file1[ num_regex.span()[1]])  
    print ('reaction number is ' + str(rxn_num))
    insertions1 = insertion_site_freq(final_sequences,template,rxn_num)
    insert_dict1 = insertions1[0] #avoiding using similar names in the workspace
    coverage = insertions1[1]


    real_insertion_list = list(insert_dict1.keys())

    add_this = 0

    if rxn_num in window1:
        add_this = base_add[0]
    elif rxn_num in window2:
        add_this = base_add[1]
    elif rxn_num in window3:
        add_this = base_add[2]
    elif rxn_num in window4:
        add_this = base_add[3]
    else:
        print('something went wrong with rxn num')
        sys.exit()
    real_insertions = [s+add_this for s in real_insertion_list] #should be able to index correctly now

    print(str(coverage)+"% coverage","total insertions: "+str(len(list(insert_dict1.keys()))))

    # fig1 = plt.figure(figsize = (30,20))
    # ax = fig1.add_subplot(1,1,1)
    # ax.scatter(real_insertions,list(insert_dict1.values()))
    # max_frequency = max(list(insert_dict1.values()))
    # figplot_scatter(ax,template,max_frequency)
    # # Append filename below as desired. 
    # fig1.savefig(output_file_prefix+'_figure.pdf')
    # should result in rxn1_828_829_F_figure.pdf as output

    # This code below is commented out because of the file name. Change as desired
    # fig1.savefig('filename.pdf')
    #
    ## Write this to a .csv file, need to write into columns instead of possible
    outp_file_loc = './CSV_Results/'
    with open(outp_file_loc+outp_name+'_results.csv','w') as file:
        # should result in rxn1_828_829_F_results.csv as output
        writer = csv.writer(file)
        writer.writerow(["insertion","count"])
        writer.writerows(zip(real_insertions,list(insert_dict1.values())))
        # should result in rxn1_828_829_F_results.csv as output
        # file.write(str(real_insertions))
        # file.write('\n')
        # file.write(str(insert_dict1.values()))
        # file.write('\n')
        # file.write('Coverage ='+str(coverage)+'%')
        file.close()

if __name__ == "__main__":
   main(sys.argv[1:])