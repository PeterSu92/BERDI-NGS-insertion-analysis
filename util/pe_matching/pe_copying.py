from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Emboss.Applications import NeedleCommandline
import re
import numpy as np
import sys

def get_coords(s):
            return ':'.join(s.description.split(' ')[0].split(':')[3:])
            
def get_sense(s):
            return s.description.split(' ')[1].split(':')[0]

def get_copied_seq(s,f_res):
            # Finds part of the template sequence by basically copying the read sequence inbetween barcode and common
            # sequence
            #return s[f_res[0].search(str(s.seq)).end():list(f_res[1].finditer(str(s.seq)))[-1].start()]
            return s[f_res[0].search(str(s.seq)).start():list(f_res[2].finditer(str(s.seq)))[-1].start()]            
            # Indices will vary!!! Must find the sequence in between the f_res sequences

def gen_copied_seq_function(f_res):
            #This does something similar to get_copied_seq but returns a function of the SeqRecord object
            # instead of a list. This is for speed then?
            return lambda s: get_copied_seq(s, f_res)
            # lambda s is just like f(s)

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
    bad_quality_reads = 0 #number of paired-end reads whose region inbetween alignment and scar has too low of a quality score to pass
    # bad_quality_reads_later = 0 #number of reads that fail the quality test after appending
    copied_too_short = 0 # if for some reason the copied region is less than the length of the smallest primer
    attempt_append = 0
    mismatched_len = 0 #number of reads that search the read and its complement badly
    read_len_list = [] #list of read lengths regardless of whether or not they pass the alignment score filter
    co_ct = 0 #number of sequences with coordinate matches
    aln_ct = 0 #number of sequences with paired end sequence matches
    append_ct = 0 #number of sequences that got appended
    #get coordinate list in the paired end reads
    count_list = []
    pe_coordL = [get_coords(s) for s in pe_seqs]
    pe_dict = {p.description:p for p in pe_seqs}
    print('begin f_seqs loop:', len(f_seqs))

    si = 0

    for s in f_seqs:
            if pe_coordL.count(get_coords(s)):
                #Apparently the above line returns a boolean so long as the count
                #isn't zero, so if the paired-end coordinates were found, the block below will be run
                co_ct += 1 
                p_index = pe_coordL.index(get_coords(s))
                pe_read = pe_seqs[p_index].reverse_complement()
                #next set of conditionals ensures the forward and paired-end reads are the same length and shaves off bases accordingly if necessary
                if len(s) > len(pe_read):
                    s = s[0:len(pe_read)]
                elif len(s) < len(pe_read):
                    pe_seqs[p_index] = pe_seqs[p_index][0:len(s)]  
                    pe_read = pe_seqs[p_index].reverse_complement()     
                if filt_seq in str(s.seq): #if the scar is present in the forward read, proceed as with the perfect match
                    copied = copied_func(s)
                    if len(copied) < 19: #since this includes a primer, should never be shorter than the shortest primer
                        copied_too_short +=1
                        continue
                    with open('temp_seq_PE.fa','w') as sh: #create temporary seq file, hopefully re-written each time
                        SeqIO.write(copied,sh,'fastq')  
                    with open('temp_temp_PE.fa','w') as PE_seq_file:
                        SeqIO.write(pe_read,PE_seq_file,'fasta')

                    needle_cline = NeedleCommandline(cmd='needle',asequence='temp_seq_PE.fa', bsequence='temp_temp_PE.fa', gapopen=10,
                                                     gapextend=0.5, outfile='PE_copied.needle') #hopefully only one needle file gets made
                    needle_cline()
                    aln_data = list(AlignIO.parse(open('PE_copied.needle'),"emboss"))
                    bin_scores = [[46,251],[213,501],[458,751],[703,1001],[952,1251],[1128,1500],[1400,1750],[1650,2000],[1800,2250],[2150,2500]] #same bin cutoff scores as alignment
                    #initialize cutoff scores
                    lo_cutoff = 0
                    hi_cutoff = 0
                    #regex searches
                    scores = score_cutoff_by_length(str(aln_data[0][0].seq).strip('-'),bin_scores)
                    lo_cutoff = scores[0]
                    hi_cutoff = scores[1]
                    match_coord_start = 0
                    match_coord_end = 0 
                    missing_align = 0
                    nonphys_overlap = 0
                    f = 0
                    if (aln_data[0].annotations['score'] >= lo_cutoff) and (aln_data[0].annotations['score'] <= hi_cutoff):
                        matched_seq_list.append(copied)
                        aln_ct += 1
                    else:
                        continue
                else: #if the scar isn't present in the forward read, start an appending process
                    with open('temp_seq_PE.fa','w') as sh: #create temporary seq file, hopefully re-written each time
                    #temp_f_seq = copied
                        SeqIO.write(s,sh,'fastq')  
                    with open('temp_temp_PE.fa','w') as PE_seq_file:
                        SeqIO.write(pe_read,PE_seq_file,'fasta')

                    needle_cline = NeedleCommandline(asequence='temp_seq_PE.fa', bsequence='temp_temp_PE.fa', gapopen=10,
                                                     gapextend=0.5, outfile='PE.needle') #hopefully only one needle file gets made
                    needle_cline()
                    aln_data = list(AlignIO.parse(open('PE.needle'),"emboss"))
                    bin_scores = [[46,251],[213,501],[458,751],[703,1001],[952,1251],[1128,1500],[1400,1750],[1650,2000],[1800,2250],[2150,2500]] #same bin cutoff scores as alignment
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
                    match_coord_start = max([bar.span()[0],bar_f.span()[0]]) #coordinates start from the first base of the forward read that aligned with the paired end read, sometimes the paired end read has bases earlier
                    match_coord_end = min([(len(aln_data[0][0].seq)-bar_f_r.span()[0]),(len(aln_data[0][0].seq)-bar_pe_r.span()[0])])# coordinates end at the last aligned base in the forward read such that no bases present on the PE read but not fwd make it
                    search_oligo = str(aln_data[0][0].seq)[match_coord_start:match_coord_end] #coordinates are currently still based off the alignment alone; search oligo is only on forward base now
                    scores = score_cutoff_by_length(search_oligo,bin_scores)
                    lo_cutoff = scores[0]
                    hi_cutoff = scores[1]
                    match_coord_start = 0
                    match_coord_end = 0 
                    missing_align = 0
                    nonphys_overlap = 0
                    f = 0
                    if (aln_data[0].annotations['score'] >= lo_cutoff) and (aln_data[0].annotations['score'] <= hi_cutoff):
                        aln_ct += 1
                    else:
                        continue
                    if len(search_oligo) < 12: #arbitrary cutoff for this aligned region based on fact that a 12 bp sequence is the smallest unique sequence in MBP
                        too_small_chunk += 1
                        continue
                    # elif len(search_oligo) > 50: #sometimes the entire region aligns, so I truncate it to just 50 bases for higher chance of alignment in the event of a mismatch surviving score filtering
                    #     search_oligo = search_oligo[len(search_oligo)-50:]
                    # #print('search oligo is '+str(len(search_oligo))+' bases long')
                    # bar1 = re.search(search_oligo,str(s.seq)) #find the aligned region in the forward sequence
                    # bar3  = re.search(str(Seq(search_oligo).reverse_complement()),str(pe_seqs[p_index].seq)) #find the aligned region's reverse complement in the actual PE sequence
                    bar4 = re.search(str(Seq(filt_seq).reverse_complement()),str(pe_seqs[p_index].seq)) # find the filt sequence's reverse complement (in this case the scar) in the actual PE 
                    # if str(type(bar3)) == "<type 'NoneType'>" or str(type(bar1)) == "<type 'NoneType'>" : #in the event there was a mismatch in the search oligo, the regex search will fail. Skip this iteration for the time being
                    #     missing_align += 1
                    #     continue
                    if str(type(bar4)) == "<type 'NoneType'>":
                        missing_pe_filt_seq += 1
                        continue
                    # elif bar4.span()[1] > bar3.span()[0]: # if some alignment happens such that part of the transposon scar aligns, this is messy and not worth dealing with
                    #     nonphys_overlap += 1
                    #     continue
                    f = quality_filter_single(pe_seqs[p_index][bar4.span()[1]:bar_pe_r],q_cutoff=20)
                    if f > 0: #if the quality of bases between the end of the aligned region and the start of the scar is good#
                        attempt_append += 1
                        temp_phred = pe_seqs[p_index].letter_annotations.values()[0][bar4.span()[1]:bar_pe_r] #temporarily dump Phred quality scores into a list
                        pe_seqs[p_index].letter_annotations = {} #clear the letter annotations so that the sequence can be changed
                        pe_seqs[p_index].seq = pe_seqs[p_index].seq[bar4.span()[1]:bar_pe_r].reverse_complement() #return only the part of the paired-end read up from end of the scar to the end of the aligned region; reverse complement
                        pe_seqs[p_index].letter_annotations = {'phred_quality':temp_phred} #now put back the new phred quality score list. Doesn't matter that it's flipped because we are only looking if any socre is below 20
                        matched_seq_list.append(pe_seqs[p_index]) #now this is appending the reverse complement of the paired end read so it can still reverse search
                        append_ct += 1 
                    else:
                        bad_quality_reads+=1
                        continue


            print si, " ", format(si/float(len(f_seqs))*100.0, '.2f'),"% percent complete            \r",
            sys.stdout.flush()
            si = si + 1
 
    # read_len_list = [f_list,pe_list]
    print ("\n\done.")
 
    # read_len_list = [f_list,pe_list]
    print ("")
    
    count_list.extend([co_ct,aln_ct]) #keep track of number of seqs with coord and align matches\
    with open('PE_statistics.csv','w') as file1:
        # should result in rxn1_828_829_F_results.csv as output
        file1.write(str(missing_filt_seq)+' reads were missing the reverse primer')
        file1.write(str(copied_too_short)+ ' reads had too small of a copied region')
        file1.write(str(missing_align)+ ' reads did not have a perfect aligned region, probably a mismatch')
        file1.write(str(too_small_chunk)+ ' reads had too small of an aligned region')
        # file1.write(str(full_align)+ ' forward reads matched all the way to the last base')
        file1.close()
    return matched_seq_list,read_len_list

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
    elif (len(sequence.strip('-')) >= 300) and (len(sequence.strip('-')) < 350):
        lo_cutoff = bin_scores[6][0]
        hi_cutoff = bin_scores[6][1]
    elif (len(sequence.strip('-')) >= 350) and (len(sequence.strip('-')) < 400):
        lo_cutoff = bin_scores[7][0]
        hi_cutoff = bin_scores[7][1]
    elif (len(sequence.strip('-')) >= 400) and (len(sequence.strip('-')) < 450):
        lo_cutoff = bin_scores[8][0]
        hi_cutoff = bin_scores[8][1]
    elif (len(sequence.strip('-')) >= 450) and (len(sequence.strip('-')) < 500):
        lo_cutoff = bin_scores[9][0]
        hi_cutoff = bin_scores[9][1]
    else:
        print(str(len(sequence.lstrip('-').strip('-')))+' is the length of the problematic read')
        raise ValueError('Sequence is either too long or too short; it is '+str(len(sequence.lstrip('-').strip('-')))+ ' bases long')
    cutoff_scores = [lo_cutoff,hi_cutoff]
    return cutoff_scores

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