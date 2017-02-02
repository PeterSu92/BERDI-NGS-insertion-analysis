# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 10:56:29 2016
Attempt to put all of Ted's ammended code in order of when I would run each step
@author: Peter Su
"""

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


# First need a function to load in fastq files
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
        seq_f = gzip.open(fpath,'rt')
        # The above line opens a gzip-compressed file in binary or text mode
        # the 'rt' indicates it is opened in text mode
    elif fpath.endswith('.fastq'):
        seq_f = open(fpath,'rt')
        # If it finds any fastq files, it simply opens them in text format
    else:
        raise ValueError('File does not end in .gz or .fastq; confirm file type')
    f_iter = SeqIO.parse(seq_f,ftype)
    return f_iter
    
    #We also need a regex function (efficient ways of searching for string) for compiling 
    #This is essentially faster since we will compile the sequences first
    
def compile_res(seqs):
            """
            Compile regex for each string in a list, return list of regex objects
            
            input:
                seqs - list of sequences you wish to filter for, e.g. CS2, transposon scar
            """
            # Takes a list of sequences you want to filter for
            # Outputs a list of regex objects that you can iterate over
            # the list comprehension below, s stands for sequences in the seqs list
            return [re.compile(s) for s in seqs]
            
#Next we want to filter through all of the fastq file to only get the sequences
#that have the corresponding required sequences in f_res
            
def filter_seqs(seqs,q_re):
    """
    Filter an iterator based on whether items match a regex object, in this case the
    sequence itself.
    inputs:
        seqs: generator object
        q_re: a regex object generated by re.compile()
    outputs:
            list of Seq objects that have sequence in them.
    """

    out_l = [s for s in seqs if q_re.search(str(s.seq))]        
            #so here we compile a list of sequences from the fastq file, stripping away
            #the other data            
            
            #List comprehension for each line in seqs??
            #text_logger.info('Finished regex filter. Kept %i sequences.',len(out_l))
    return out_l
            
            #Next we will start filtering based on the sequences desired in compile_res
            
def filter_sample(f_name,pe_name,template,f_filt_seqs,r_filt_seqs):
        """
        Outputs filtered sequences as "dictionary", indexed by barcode        
        Sequences will be aligned to the provided template.
        Parts of the template not represented will be '-'
        
        inputs:
            fname - str; file name
            pe_name - str, paired-end read file name
            templates - str, my guess is sequence of template to align to
            f_filt_seq - list of sequences to filter for
            r_filt_seq - list of sequences to filter for on the paired-end side
        """
        
        # Compile regexes
        f_res = compile_res(f_filt_seqs)
        pe_res = compile_res(r_filt_seqs)
        
        # Load FASTQ files as generator
        f_seqs = load_ngs_file(f_name)
        pe_seqs = load_ngs_file(pe_name)
        # Create lists of SeqRecord objects from the generators
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
        #The outermost list will be the size of f_filt_seqs (same as r_filt_seqs)
        #The inner lists ought to contain the SeqRecord reads that had that particular
        #sequence we were searching for (i.e. CS or TR). The f_seqs2 list will contain one of the CS's,
        #while the pe_seqs2 will contain the other (so indices will be messed up until I figure out which
        #read consistently has which CS). TLDR take indices in this next part with a grain of salt, and it's
        #overall kind of clunky. 

        f_seqs3 = f_seqs2[2]

        # Repeat for paired-end reads

        pe_seqs3 = pe_seqs2[1]

        print(str(len(f_seqs3))+' forward reads have the sequence of interest (MBP forward primer)')
        print(str(len(pe_seqs3))+' paired-end reads have the sequence of interest (transposon scar)')
        # Now that only sequences containing BOTH the CS and the TR have been filtered for,
        # the paired-end matching can occur
        
        seqs = filter_pe_mismatch(f_seqs3,pe_seqs3,gen_copied_seq_function(f_res),f_filt_seqs[1]) #forward reads use the MBP reverse primer here
        with open('matched_seq_PE.fa','w') as sh: #create temporary seq file, hopefully re-written each time
                #temp_f_seq = copied
             SeqIO.write(seqs,sh,'fastq')
        # read_len_postalign_list = s1[1]
        print(str(len(seqs))+' forward reads have a paired-end match')
        seqs = quality_filter(seqs,q_cutoff=20)
        print(str(len(seqs))+' forward reads survived the Phred score quality filter')
        # Note this returns only the forward sequences
        # Now the final step of aligning and filtering by alignment scores happens
        #run the alignment to the template file, varying cutoff by length 
        #initialize bins
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
        
def quality_filter(seqs,q_cutoff=20):
            """
            removes sequences that have any bases below the cutoff score specified
            
            inputs:
                seqs
                q_cutoff - int, score below which a sequence will be removed
                
            outputs:
                out_l = list of sequences that survived the filtering
            """
            out_l = [s for s in seqs
                    if not any(s.letter_annotations['phred_quality'] < np.ones(len(s.letter_annotations['phred_quality']))*q_cutoff)]
            return out_l

def quality_filter_single(seqs,q_cutoff=20):
            """
            removes sequences that have any bases below the cutoff score specified
            
            inputs:
                seqs
                q_cutoff - int, score below which a sequence will be removed
                
            outputs:
                out_l = list of sequences that survived the filtering
            """
            out_l = 0
            if not any(seqs.letter_annotations['phred_quality'] < np.ones(len(seqs.letter_annotations['phred_quality']))*q_cutoff):
                out_l += 1
            return out_l
            
    # This function throws out any aligned reads below a certain cutoff score
def cull_alignments(aln_data, lo_cutoff, hi_cutoff):
        """
        Filters out reads that returned too low of a score from the needle command
        inputs:
                aln_data - AlignIO parse/read-in of the alignments
                lo_cutoff - int, low score cutoff based on needle
                hi_cutoff, int, high score cutoff based on needle
        outputs:
            new_seqs - list of sequences that survive the cutoff
        """
        new_seqs = []
        for alignment in aln_data:
            if (alignment.annotations['score'] >= lo_cutoff) and (alignment.annotations['score'] <= hi_cutoff):
                        #Template should have no gaps and should contain the whole
                        # non-template sequence
                            joined_align = [r for r,t in zip(alignment[1],alignment[0]) if t != '-']
                            new_read = SeqRecord(''.join(joined_align))
                            new_seqs.append(new_read)
                            new_seqs[-1].annotations['alnscore'] = alignment.annotations['score']
        return new_seqs
        
def alignment_filter(seqs,template, lo_cutoff,hi_cutoff, bin_num, gapopen = 10, gapextend = 0.5):
    '''
         Aligns sequences to a template file and filters the aligned reads by score
         
         Inputs:
             seqs - list of SeqRecord objects
             template - str, DNA sequence of the template of interest
         Outputs:
             new_seqs - list of SeqRecord objects that survived the filtering
    '''
    # Save the template and sequences as temporary fasta files 
    template_f_name = 'temptemplate.fa'
    seqs_f_name = 'tempseq.fa'
            
    with open(seqs_f_name,'w') as sh:
        SeqIO.write(seqs,sh,'fastq')
                
    with open(template_f_name,'w') as temp_seq_file:
        # make temp sequence file for alignment
        temp_seq = SeqRecord(Seq(template),id='template',name = 'template')
        SeqIO.write(temp_seq,temp_seq_file,'fasta')
                
    # Generate alignment command, run the alignment
    ofilen = bin_num+'temp_'+str(uuid.uuid4())+'.needle' #needle files have bin number
    needle_cline = NeedleCommandline(asequence=template_f_name, bsequence=seqs_f_name, gapopen=gapopen,
                                             gapextend=gapextend, outfile=ofilen)
    needle_cline()
            
    aln_data = AlignIO.parse(open(ofilen),"emboss")
    # for some reason this generator stuff is not working for me so I just set it as a list
    aln_data_list = list(aln_data)
    new_seqs = cull_alignments(aln_data_list, lo_cutoff, hi_cutoff)
            
    return new_seqs
        
def get_coords(s):
            return ':'.join(s.description.split(' ')[0].split(':')[3:])
            
def get_sense(s):
            return s.description.split(' ')[1].split(':')[0]
        
def get_copied_seq(s,f_res):
            # Finds part of the template sequence by basically copying the read sequence inbetween barcode and common
            # sequence
            #return s[f_res[0].search(str(s.seq)).end():list(f_res[1].finditer(str(s.seq)))[-1].start()]
            return s[f_res[2].search(str(s.seq)).end():list(f_res[1].finditer(str(s.seq)))[-1].end()]            
            # Indices will vary! Must find the sequence in between the f_res sequences
            # Don't understand finditer super well
# def get_copied_seq_half(s,f_res):
#             # Finds part of the template sequence by basically copying the read sequence inbetween barcode and common
#             # sequence
#             #return s[f_res[0].search(str(s.seq)).end():list(f_res[1].finditer(str(s.seq)))[-1].start()]
#             return s[f_res[0].search(str(s.seq)).start():]  
def gen_copied_seq_function(f_res):
            #This does something similar to get_copied_seq but returns a function of the SeqRecord object
            # instead of a list. This is for speed then?
            return lambda s: get_copied_seq(s, f_res)
            # lambda s is just like f(s)

def score_cutoff_by_length(sequence,bin_scores):

    """
    Determines the size of the sequence and assigns a corresponding bin low and high cutoff score for the Needleman-Wunsch algorithm

    Inputs:
        sequence - str, the characters (DNA bases here) for which a score cutoff is determined
        bin_scores - list of lists of integers containing the corresponding scores

    Outputs:
        cutoff_scores - list of two integers, the low and high cutoff scores

    """
    #initialize variables
    lo_cutoff = 0
    hi_cutoff = 0
    # determine length and set scores
    if len(sequence.lstrip('-').strip('-')) < 50 :
        lo_cutoff = bin_scores[0][0]
        hi_cutoff = bin_scores[0][1]
    elif len(sequence.lstrip('-').strip('-')) >= 50 and len(sequence.lstrip('-').strip('-')) < 100:
        lo_cutoff = bin_scores[1][0]
        hi_cutoff = bin_scores[1][1]
    elif len(sequence.lstrip('-').strip('-')) >= 100 and len(sequence.lstrip('-').strip('-')) < 150:
        lo_cutoff = bin_scores[2][0]
        hi_cutoff = bin_scores[2][1]
    elif len(sequence.lstrip('-').strip('-')) >= 150 and len(sequence.lstrip('-').strip('-')) < 200:
        lo_cutoff = bin_scores[3][0]
        hi_cutoff = bin_scores[3][1]
    elif len(sequence.lstrip('-').strip('-')) >= 200 and len(sequence.lstrip('-').strip('-')) < 250:
        lo_cutoff = bin_scores[4][0]
        hi_cutoff = bin_scores[4][1]
    elif len(sequence.lstrip('-').strip('-')) >= 250 and len(sequence.lstrip('-').strip('-')) <= 301:
        lo_cutoff = bin_scores[5][0]
        hi_cutoff = bin_scores[5][1]
    else:
        print(str(len(sequence.lstrip('-').strip('-')))+' is the length of the problematic read')
        raise ValueError('Sequence is either too long or too short; it is '+str(len(sequence.lstrip('-').strip('-')))+ ' bases long')
    cutoff_scores = [lo_cutoff,hi_cutoff]
    return cutoff_scores
            
def filter_pe_mismatch(f_seqs,pe_seqs,copied_func,filt_seq): #Now edited to use the Needleman-Wunsch algorithm for paired-end filtering.
    """
    
    Inputs:
        f_seqs - list of sequences that survived filtering so far (SeqRecord objects)      
        pe_seqs - list of paired end sequences that survived filtering so far        
        copied_func - list of sequences; output of gen_copied_seq_function(f_res)
        
    Output:
        matched_seq_list - list of sequences that had paired end matches
    """
    #initialize variables
    matched_seq_list = []
    f_list = []
    pe_list = []
    missing_filt_seq = 0 #number of forward reads missing the back filter sequence
    missing_pe_filt_seq = 0 #number of PE reads missing the scar (shouldn't happen?)
    too_small_chunk = 0 #number of forward reads that had too small of a chunk to be kept
    full_align = 0 #number of reads where the last base of the forward read aligned to the PE read
    copied_too_short = 0 # if for some reason the copied region is less than the length of the smallest primer
    read_len_list = [] #list of read lengths regardless of whether or not they pass the alignment score filter
    co_ct = 0 #number of sequences with coordinate matches
    aln_ct = 0 #number of sequences with paired end sequence matches
    append_ct = 0 #number of sequences that got appended
    count_list = []
       #get coordinate list in the paired end reads
    pe_coordL = [get_coords(s) for s in pe_seqs]
    pe_dict = {p.description:p for p in pe_seqs}
    #pe_coord_dict = {pe_coordL[s]:pe_dict}
    print('begin f_seqs loop:', len(f_seqs))

    si = 0

    for s in f_seqs:
            if pe_coordL.count(get_coords(s)):
                #Apparently the above line returns a boolean so long as the count
                #isn't zero, so if the paired-end coordinates were found, the block below
                # will be run
                co_ct += 1  
                p_index = pe_coordL.index(get_coords(s))
                pe_read = pe_seqs[p_index].reverse_complement()
                if len(s) > len(pe_read):
                    s = s[0:len(pe_read)]
                elif len(s) < len(pe_read):
                    pe_seqs[p_index] = pe_seqs[p_index][0:len(s)]  
                    pe_read = pe_seqs[p_index].reverse_complement()     

                if filt_seq in str(s.seq): #if the scar is present in the forward read, proceed as with the perfect match
                    #Get part of the sequence that was actually copied
                    copied = copied_func(s)
                    if len(copied) < 19:
                        copied_too_short +=1
                        continue
                    with open('temp_seq_PE.fa','w') as sh: #create temporary seq file, hopefully re-written each time
                        SeqIO.write(copied,sh,'fastq')  
                    with open('temp_temp_PE.fa','w') as PE_seq_file:
                        SeqIO.write(pe_read,PE_seq_file,'fasta')

                    needle_cline = NeedleCommandline(asequence='temp_seq_PE.fa', bsequence='temp_temp_PE.fa', gapopen=10,
                                                     gapextend=0.5, outfile='PE.needle') #hopefully only one needle file gets made
                    needle_cline()
                    aln_data = list(AlignIO.parse(open('PE.needle'),"emboss"))
                    bin_scores = [[46,251],[213,501],[458,751],[703,1001],[952,1251],[1128,1500]] #same bin cutoff scores as alignment
                    #initialize cutoff scores
                    lo_cutoff = 0
                    hi_cutoff = 1500
                    cond_bar = re.search('[AGCT]+',str(aln_data[0][1].seq)) #search forwards through reverse complement of PE read, find first base that aligned.
                    cond_bar_f = re.search('[AGCT]+',str(aln_data[0][0].seq)[-1:0:-1]) #search backwards through forward read, find first base that aligned
                    cond_match_coord_start =cond_bar.span()[0] #coordinates start from the first base of the forward read that aligned with the paired end read
                    cond_match_coord_end = len(aln_data[0][0].seq)-cond_bar_f.span()[0]
                    cond_search_oligo = str(aln_data[0][0].seq)[cond_match_coord_start:cond_match_coord_end]
                    scores = score_cutoff_by_length(cond_search_oligo,bin_scores)
                    lo_cutoff = scores[0]
                    hi_cutoff = scores[1]
                    match_coord_start = 0
                    match_coord_end = 0 
                    match_len = 0
                    missing_align = 0
                    nonphys_overlap = 0
                    f = 0
                    if (aln_data[0].annotations['score'] >= lo_cutoff) and (aln_data[0].annotations['score'] <= hi_cutoff):
                        matched_seq_list.append(copied)
                        aln_ct += 1
                    else:
                        continue
                else:
                    with open('temp_seq_PE.fa','w') as sh: 
                        SeqIO.write(s,sh,'fastq')  
                    with open('temp_temp_PE.fa','w') as PE_seq_file:
                        SeqIO.write(pe_read,PE_seq_file,'fasta')

                    needle_cline = NeedleCommandline(asequence='temp_seq_PE.fa', bsequence='temp_temp_PE.fa', gapopen=10,
                                                     gapextend=0.5, outfile='PE.needle') #hopefully only one needle file gets made
                    needle_cline()
                    aln_data = list(AlignIO.parse(open('PE.needle'),"emboss"))
                    bin_scores = [[46,251],[213,501],[458,751],[703,1001],[952,1251],[1128,1500]] #same bin cutoff scores as alignment
                    #initialize cutoff scores
                    lo_cutoff = 0
                    hi_cutoff = 1500
                    #if the scar isn't found on the forward read
                    missing_filt_seq +=1
                    bar = re.search('[AGCT]+',str(aln_data[0][1].seq)) #search forwards through reverse complement of PE read, find first base that aligned.
                    bar_pe_r = re.search('[AGCT]+',str(aln_data[0][1].seq)[-1:0:-1]) #search backwards through PE read, find first base that aligned
                    bar_f = re.search('[AGCT]+',str(aln_data[0][0].seq))
                    bar_f_r = re.search('[AGCT]+',str(aln_data[0][0].seq)[-1:0:-1]) # search backwards through fwd read to find the first base that aligned
                    # if bar.span()[0] < bar_f.span()[0]: #this is if the fwd strand search goes longer until it hits an aligned base
                    match_coord_start = max(bar.span()[0],bar_f.span()[0]) #coordinates start from the first base of the forward read that aligned with the paired end read, sometimes the paired end read has bases earlier
                    match_coord_end = min(len(aln_data[0][0].seq)-bar_f_r.span()[0],len(aln_data[0][0].seq-bar_pe_r.span()[0])# coordinates end at the last aligned base in the forward read such that no bases present on the PE read but not fwd make it
                    search_oligo = str(aln_data[0][0].seq)[match_coord_start:match_coord_end] #coordinates are currently still based off the alignment alone; search oligo is only on forward base now
                    scores = score_cutoff_by_length(search_oligo,bin_scores)
                    lo_cutoff = scores[0]
                    hi_cutoff = scores[1]
                    if (aln_data[0].annotations['score'] >= lo_cutoff) and (aln_data[0].annotations['score'] <= hi_cutoff):
                        aln_ct += 1
                    else:
                        continue
                    if len(search_oligo) < 12: # this number is arbitrary right now
                        too_small_chunk += 1
                        continue
                    bar1 = re.search(search_oligo,str(s.seq)) #find the aligned region in the forward sequence
                    # f = quality_filter_single(pe_seqs[p_index][bar4.span()[0]:bar3.span()[0]],q_cutoff=20)
                    # if f > 0: #if the quality of bases between the end of the aligned region and the start of the scar is good#
                    #     bar2 = pe_read_rev.find(filt_seq) # find the filter sequence in the reverse complement of the PE read, for the purpose of appending a region
                    if str(type(bar1)) == "<type 'NoneType'>": #if the aligned region can't be found in the forward read
                        missing_align += 1
                        continue
                    elif bar1.span()[1] > len(s.seq):
                         raise ValueError('End coordinate is '+str(bar1.span()[1])+' but sequence is '+str(len(s.seq)))
                    elif bar1.span()[1] == len(s.seq):
                        full_align += 1
                        # print ("bar5.span()[1] is "+str(bar5.span()[1])+ " bar2.span()[0] is "+ str(bar2.span()[0]))
                        # pe_append = pe_read_rev[bar5.span()[1]:bar2] #hopefully this returns the part of the paired-end read from the last base of alignment to the scar
                        # attempt_append += 1
                    temp_phred = s.letter_annotations.values()[0][0:bar1.span()[1]] #temporarily dump Phred quality scores into a list
                    s.letter_annotations = {} #clear the letter annotations so that the sequence can be changed
                    s.seq = s.seq[0:bar1.span()[1]] #return only the part of the forward read up to the end of the aligned region. This way, no junk gets kept in the case of a short read
                    s.letter_annotations = {'phred_quality':temp_phred} #now put back the new phred quality score list
                    matched_seq_list.append(s)

            print si, " ", format(si/float(len(f_seqs))*100.0, '.2f'),"% percent complete            \r",
            sys.stdout.flush()
            si = si + 1
    print ("done.")
    print ("")
    
    count_list.extend([co_ct,aln_ct]) #keep track of number of seqs with coord and align matches
    # print(str(co_ct)+' forward reads had the coordinates of the PE read nearby')
    print(str(missing_filt_seq)+' reads were missing the reverse primer')
    print(str(copied_too_short)+ ' reads had too small of a copied region')
    # print(str(missing_pe_filt_seq)+ 'paired-end reads are missing the scar!!!!')
    print(str(missing_align)+ ' reads did not have a perfect aligned region')
    print(str(too_small_chunk)+ ' reads had too small of an aligned region')
    print(str(full_align)+ ' forward reads matched all the way to the last base')
    # print(str(nonphys_overlap)+' reads had part of the scar or some nonphysical overlap')
    # print(str(bad_quality_reads_first)+ 'reads had poor quality in the region to be appended')
    # print(str(bad_quality_reads_later)+' reads had overall poor quality in the final sequence')
    # print(str(attempt_append)+ ' reads started to be appended before perhaps failing in the added bases quality')
    # print(str(append_ct)+ ' reads had appended parts from the paired-end read')
    return matched_seq_list
    
#Insertion site functions and code
def insertion_chunks(final_seqs):
    '''
    Should create a list of contiguous stretches of DNA for each sequence in
    new_seqs
    
    Input:
        final_seqs - list of SeqRecord objects containing the filtered sequences from the 
        EMBOSS needle alignment algorithm. 
        template - str, DNA sequence of the template
        
    Outputs:
        chunk_dict - dict of chunks of continuous DNA segments per read.
        Each entry has key 'insert_site' with a corresponding value of a list of
        the lengths of each continuous chunk (starting from the first aligned bit)
        insertions - list of positions where the first chunk of continuous DNA began
        for each particular read. Will have redundant entries in most cases, as
        our method often results in multiple insertions at any given site.
    '''
    chunk_dict = {}
    insertions = []
    chunk_size = 0
    max_chunks = 0
    reads_at_end = 0
    large_chunk_reads = 0
    end_dashes = 0
    perfect_matches = 0
    max_chunks_exceeded = 0
    other_scenario = 0
    discarded_reads = 0

    for i in range(len(final_seqs)):
           end_pos = 0 #forward search starts at the beginning
           #insert_site = 0
           num_chunks = 0
           seq_chunks = []
           insert_site = 0
           
           total_len = 0
           #print('Current sequence: ' +str(i+1)) #keep this only for test sequences
           if str(final_seqs[i].seq)[-1] == '-':
              discarded_reads += 1
              end_dashes +=1
              continue
           while total_len < len(final_seqs[i].seq):
              bar=re.search('[AGCT]+',str(final_seqs[i].seq)[end_pos:-1:1]) #forward search: from start to finish
              if str(type(bar)) == "<type 'NoneType'>":
                   #If this happens, we'll know the end was reached without finding a suitable insertion
                    reads_at_end += 1
                    break

              if end_pos >= 300: #this prevents a nonphysical insertion from happening
                    end_pos = end_pos-4
                    insertions.append(insert_site)
                    break
              elif abs((bar.span()[1]-bar.span()[0])+1) == (len(final_seqs[i].seq.lstrip('-').strip('-'))): #perfect match occurs
                     insert_site = bar.span()[0] #forward search stops at the first base of the DNA chunk
                     insertions.append(insert_site)
                     chunk_dict.update({insert_site:seq_chunks})
                     seq_chunks.append(bar.span()[1]-bar.span()[0])
                     perfect_matches +=1
                     break
              elif abs((bar.span()[1]-bar.span()[0])) <= chunk_size: #if a chunk is small enough, set index correspondingly but keep searching through the alignment
                     num_chunks += 1
                     end_pos += bar.span()[1]
                     insert_site += bar.span()[1]
                     span_length = abs(bar.span()[1]-bar.span()[0])
                     seq_chunks.append(span_length)
                     total_len += bar.span()[1]
                     continue
              elif (abs((bar.span()[1]-bar.span()[0])) > chunk_size) and (abs((bar.span()[1]-bar.span()[0])) != len(final_seqs[i].seq.lstrip('-'))-total_len): #if chunk too large, get rid of the alignment
                     discarded_reads += 1
                     large_chunk_reads +=1
                     break
              elif len(final_seqs[i].seq.strip('-')) != len(final_seqs[i].seq.lstrip('-')): #gets rid of alignments with gaps at the 3' end
                     discarded_reads +=1
                     break
              elif num_chunks > max_chunks: #too many chunks leads to an alignment being thrown out. 
                     discarded_reads += 1
                     max_chunks_exceeded +=1
                     break
         
              else: #This should not happen now, but if it does, will document it
                    other_scenario +=1
                    break

    print (str(reads_at_end)+ ' reads reached the end without a suitable insertion')    
    print (str(discarded_reads)+' reads were discarded :(')
    print(str(large_chunk_reads)+' reads had too large of a chunk')
    print(str(max_chunks_exceeded)+' reads had too many chunks')
    print(str(end_dashes)+' reads had end dashes')
    print(str(perfect_matches)+' reads are perfect matches')
    print(str(other_scenario) +' reads did not satisfy any of the criteria')


    return chunk_dict, insertions
    
def insertion_site_freq(final_seqs,template,reaction_number):
    '''
    Calculates the frequency of insertions at all sites in a template sequence
    
    Inputs:
        final_seqs - list of SeqRecord objects containing the filtered sequences from the 
        EMBOSS needle alignment algorithm.
        template - str, DNA sequence of the template
    
    Ouputs:
        insert_dict: dictionary of insertion site (nucleotide):insertion counts 
            in the template sequence for each site that has at least one insertion
        coverage: int,% of all possible sites that have at least one insertion
    '''
    
    chunky_dict = insertion_chunks(final_seqs)
    insert_dict = {} #dict consisting of insertion site: insertion count pairs
    insertion_list = chunky_dict[1]
    #insertion_frequencies = []
    insertion_site_list = list(set(insertion_list))
    
    for sites in insertion_site_list:
        insert_dict[sites] = insertion_list.count(sites)
        #insertion_frequencies.append(insertion_list.count(sites))
    
    coverage = 100*len(list(insertion_site_list))/len(template)
    
    return insert_dict,coverage


def figplot_scatter(ax,template,max_frequency):
    '''
    Set of simple commands to make the scatterplot figure look nice
    
    Inputs: 
        -ax, set of individual plots
        -template - str, template sequence
        -max_frequency - maximum frequency
    
    '''
    ax.set_xlim([0,1116])
    ax.set_yscale('log')
    ax.set_ylim([0,max_frequency])
    ax.set_xlabel('Insertion site',fontsize = 30)
    ax.set_ylabel('Frequency',fontsize = 30)
    #ax.legend(fontsize=20)
    



#### Actually running the code ######
import sys, getopt, os

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

    print(str(coverage)+"% coverage","total insertions"+str(len(list(insert_dict1.keys()))))

    fig1 = plt.figure(figsize = (30,20))
    ax = fig1.add_subplot(1,1,1)
    ax.scatter(real_insertions,list(insert_dict1.values()))
    max_frequency = max(list(insert_dict1.values()))
    figplot_scatter(ax,template,max_frequency)
    # Append filename below as desired. 
    fig1.savefig(output_file_prefix+'_figure.pdf')
    # should result in rxn1_828_829_F_figure.pdf as output

    # This code below is commented out because of the file name. Change as desired
    # fig1.savefig('filename.pdf')
    #
    ## Write this to a .csv file, need to write into columns instead of possible
    outp_file_loc = '../CSV_Results/'
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
##
