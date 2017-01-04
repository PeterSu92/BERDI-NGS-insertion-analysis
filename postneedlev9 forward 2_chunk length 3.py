# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 10:41:09 2016
Subset of the original script to filter only the needle files
@author: PSu92_000
"""

import os
import re
import sys
import time
import uuid
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import csv
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Emboss.Applications import NeedleCommandline




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
            if (alignment.annotations['score'] >= lo_cutoff) and (alignment.annotations['score'] < hi_cutoff):
                        #Template should have no gaps and should contain the whole
                        # non-template sequence
                #if not str(alignment[0].seq).count('-') > 0:
                            joined_align = [r for r,t in zip(alignment[1],alignment[0]) if t != '-']
                            new_read = SeqRecord(''.join(joined_align))
                            new_seqs.append(new_read)
                            new_seqs[-1].annotations['alnscore'] = alignment.annotations['score']
        return new_seqs
        
#def alignment_filter(seqs,template, gapopen = 10, gapextend = 0.5, lo_cutoff = 300,hi_cutoff=1000):
#    '''
#         Aligns sequences to a template file and filters the aligned reads by score
#         
#         Inputs:
#             seqs - list of SeqRecord objects
#             template - str, DNA sequence of the template of interest
#         Outputs:
#             new_seqs - list of SeqRecord objects that survived the filtering
#    '''
#    # Save the template and sequences as temporary fasta files 
#    template_f_name = 'temptemplate.fa'
#    seqs_f_name = 'tempseq.fa'
#            
#    with open(seqs_f_name,'w') as sh:
#        SeqIO.write(seqs,sh,'fastq')
#                
#    with open(template_f_name,'w') as temp_seq_file:
#        # make temp sequence file for alignment
#        temp_seq = SeqRecord(Seq(template),id='template',name = 'template')
#        SeqIO.write(temp_seq,temp_seq_file,'fasta')
#                
#    # Generate alignment command, run the alignment
#    ofilen = 'temp_'+str(uuid.uuid4())+'.needle'
#    needle_cline = NeedleCommandline(asequence=template_f_name, bsequence=seqs_f_name, gapopen=gapopen,
#                                             gapextend=gapextend, outfile=ofilen)
#    needle_cline()
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
    chunk_size = 3
    max_chunks = 2
    read_at_ends = 0
# if reaction_number in reverse_search:
#     end_pos_default = -1
#     final_pos = 0
# elif reaction_number in forward_search:
#     end_pos_default = 0
#     final_pos = -1
# else:
#     print('Error your numbering is terrible')
    discarded_reads = 0

    for i in range(len(final_seqs)):
           end_pos = 0 #forward search starts at the beginning
           #insert_site = 0
           num_chunks = 0
           seq_chunks = []
           insert_site = 0
           
           total_len = 0
           print('Current sequence: ' +str(i+1)) #keep this only for test sequences
           if str(final_seqs[i].seq)[-1] == '-':
              discarded_reads += 1
              print('Sequence '+str(i+1)+ ' had dashes at the 3prime end')
              continue
           while total_len < len(final_seqs[i].seq):
              bar=re.search('[AGCT]+',str(final_seqs[i].seq)[end_pos:-1:1]) #forward search: from start to finish
              if str(type(bar)) == "<type 'NoneType'>":
                   #If this happens, we'll know the last base of the previous
                   # chunk was the insertion site, so we set it as such here.                
                print('Sequence '+str(i+1)+'end reached')
                    #insert_site = end_pos
                if end_pos >= 300: #this prevents a nonphysical insertion from happening
                    end_pos = end_pos-4
                    insertions.append(insert_site)
                break
              elif abs((bar.span()[1]-bar.span()[0])+1) == (len(final_seqs[i].seq.lstrip('-').strip('-'))): #perfect match occurs
                     insert_site = bar.span()[0] #forward search stops at the first base of the DNA chunk
                     insertions.append(insert_site)
                     chunk_dict.update({insert_site:seq_chunks})
                     seq_chunks.append(bar.span()[1]-bar.span()[0])
                     print('Sequence '+str(i+1)+' had a perfect match')
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
                     print('Sequence '+str(i+1)+' had too large of a chunk')
                     break
              elif len(final_seqs[i].seq.strip('-')) != len(final_seqs[i].seq.lstrip('-')): #gets rid of alignments with gaps at the 3' end
                     discarded_reads +=1
                     break
              elif num_chunks > max_chunks: #too many chunks leads to an alignment being thrown out. 
                     discarded_reads += 1
                     break
                    # elif abs((bar.span()[1]-bar.span()[0])) <=8: #if the contiguous chunk of DNA is too small, move on but save the indexing and the length searched so that the
                        # insert site isn't miscaculated, as the next iteration of the search starts the regex indexing back at 0
                        #print('an insertion had a small chunk;'+ ' end pos: '+str(bar.span()[1]))
                        # end_pos += bar.span()[1]
                        # #we add the entire length of the spanned region, through the chunk of DNA spanned, to the insert site
                        # insert_site += (bar.span()[1])
                        # span_length = abs(bar.span()[1]-bar.span()[0])
                        # seq_chunks.append(span_length)
                        # total_len += bar.span()[1]
                        # num_chunks +=1
                        # continue
                    # elif ((bar.span()[1]-bar.span()[0]) >=8) and (bar.span()[0] == 0): #this shouldn't ever really happen for most PCRs in this category. 
                    #     insert_site += bar.span()[0]
                    #     print('Insertion at start of template') #printing this will tell us if this happens
                    #     num_chunks += 1
                    #     insertions[i] = insert_site
                    #     span_length = abs(bar.span()[1]-bar.span()[0])
                    #     seq_chunks.append(span_length)
                    #     total_len += span_length
                    #     break
         
              else: #This should not happen now, but if it does, it sets the insert site as the FIRST base of the contiguous region, which is why bar.span()[0] is used
                    #insert_site = bar.span()[0]
                    print('Sequence '+str(i+1)+ 'had an insertion weirdly')
                    if bar.span()[0] >= 300:
                       insert_site = bar.span()[0]-4
                       insertions.append(insert_site)
                       chunk_dict.update({insert_site:seq_chunks})
                    else:
                       insert_site = bar.span()[0]
                       insertions.append(insert_site)
                       chunk_dict.update({insert_site:seq_chunks})

                    break
                #if insert_site != 0:
            
    # insertions.append(insert_site)
            

    
    
    return chunk_dict, insertions
    
def insertion_site_freq(final_seqs):
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
    
    coverage = 100*len(list(insertion_site_list))/1116 
    
    return insert_dict,coverage
    
## Part that gets run ###


#template_file = '../temptemplate.fa'
#template = str(list(SeqIO.parse(template_file,'fasta'))[0].seq)
### Append this prefix as desired ###
output_file_prefix = 'PCR3 4 chunk'
#Bin cutoff scores
bin_scores = [[42,251],[205,501],[446,751],[687,1001],[928,1251],[1100,1500]]
final_sequences = []
newseqs = []
needle_files = glob.glob('*.needle')
for i in range(len(needle_files)):
    aln_data = AlignIO.parse(open(needle_files[i]),"emboss")
    # for some reason this generator stuff is not working for me
    aln_data_list = list(aln_data)
    #print('Sequences in bin '+str(i+1)+' before alignment: '+str(len(aln_data_list)))
    new_seqs = cull_alignments(aln_data_list, lo_cutoff=bin_scores[i][0], hi_cutoff=bin_scores[i][1])
    print('Sequences in bin '+str(i+1)+' after alignment: '+str(len(new_seqs)))
    newseqs.append(new_seqs)

    insertions2 = insertion_site_freq(new_seqs)
    insert_dict2 = insertions2[0] #avoiding using similar names in the workspace
    real_insertion_list = list(insert_dict2.keys())
    # print('Insertions in bin '+str(i+1)+': '+str(len(real_insertion_list)))
    if len(real_insertion_list) == 0:
        print('Insertions in bin '+str(i+1)+': 0')
    else:
        print('Insertions in bin '+str(i+1)+': '+str(len(real_insertion_list)))
#Coverage amended here to show coverage over entire sequence of MBP
    #coverage = insertions1[1]


#merge sublits in new_seqs into one master list
final_sequences = [item for sublist in newseqs for item in sublist]
len_list = [len(str(s.seq)) for s in final_sequences]
print ('min_length ' + str(min(len_list)) + '; ' + 'max_length' + str(max(len_list)))
insertions1 = insertion_site_freq(final_sequences)
insert_dict1 = insertions1[0] #avoiding using similar names in the workspace
#Coverage amended here to show coverage over entire sequence of MBP
coverage = insertions1[1]

real_insertion_list = list(insert_dict1.keys())
#varies depending on portion of MBP used as template
real_insertions = [s+304 for s in real_insertion_list]
print(str(coverage)+"% coverage","total insertions "+str(len(list(insert_dict1.keys()))))
with open(output_file_prefix+'_results.csv','w') as file:
        # should result in rxn1_828_829_F_results.csv as output
        writer = csv.writer(file)
        writer.writerow(["insertion","count"])
        writer.writerows(zip(real_insertions,list(insert_dict1.values())))
    #if cleanup:
    #   logging.info('Cleaning up temp files')
    #   os.remove(template_fname)
    #   os.remove(seqs_fname)
    #   os.remove(ofilen)
