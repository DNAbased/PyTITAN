import pandas as pd
import numpy as np
import sys
from pyfaidx import Fasta
import os


class ConversionHandler:
    def __init__(self):
        pass


    @staticmethod
    def pfm_to_pwm_folder(pfm_paths, pwm_folder, pseudocount=0.8, bg={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}):
        '''Convert all PFMs in a folder to PWMs in a separate folder'''
        for i in pfm_paths:
            pfm = i.split('/')[-1].split('_')
            pfm_id = pfm[0]
            pfm_name = pfm[1]
            pfm_strand = pfm[2].split('.')[0]
            pfm_matrix = pd.read_csv(i, sep='\t', header=None, index_col=0)
            for i in ['A', 'C', 'G', 'T']:
                pfm_matrix.loc[i] = pfm_matrix.loc[i] + (bg[i] * pseudocount)
            pfm_matrix = pfm_matrix / pfm_matrix.sum()
            for i in ['A', 'C', 'G', 'T']:
                pfm_matrix.loc[i] = pfm_matrix.loc[i] / bg[i]
            pfm_matrix = np.log2(pfm_matrix)
            pfm_matrix.to_csv(pwm_folder + pfm_id + '_' + pfm_name + '_' + pfm_strand + '.pwm')


    @staticmethod
    def pfm_to_icm_folder(pfm_paths, info_cont_folder, pseudocount=0.8, bg={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}):
        '''Convert all PFMs in a folder to ICMs in a separate folder'''
        for i in pfm_paths:
            pfm = i.split('/')[-1].split('_')
            pfm_id = pfm[0]
            pfm_name = pfm[1]
            pfm_strand = pfm[2].split('.')[0]
            pfm_matrix = pd.read_csv(i, sep='\t', header=None, index_col=0)
            for i in ['A', 'C', 'G', 'T']:
                pfm_matrix.loc[i] = pfm_matrix.loc[i] + (bg[i] * pseudocount)
            icm_matrix = pfm_matrix / pfm_matrix.sum()
            icm_matrix_log = np.log2(icm_matrix)
            icm_matrix_sum = (icm_matrix * icm_matrix_log).sum() + 2
            icm_matrix = icm_matrix * icm_matrix_sum
            icm_matrix.to_csv(info_cont_folder + pfm_id + '.icm', header=False, index=False)


class DataHandler:
    def __init__(self):
        pass


    @staticmethod
    def get_files(filepath):
        '''Get all files and their path from a folder'''
        f = []
        for path, dirs, files in os.walk(filepath):
            for file in files:
                f.append(os.path.join(path, file))
        return(f)


    @staticmethod
    def create_id(df):
        '''Add a (unique) ID column to each variant of a data frame'''
        for i in df.index:
            df.at[i, 'ID'] = str(df.at[i, 'CHROM']) + '_' + str(df.at[i, 'POS']) + '_' + df.at[i, 'REF'] + '_' + df.at[i, 'ALT']
        return(df)


    @staticmethod
    def get_seq(df, genome='hs37d5.fa', n=50):
        '''Get the genomic sequence of a position including up-/downstream bases'''
        for i in df.index:
            df.at[i, 'SEQ'] = Fasta(genome)[str(df.at[i, 'CHROM'])][(df.at[i, 'POS'] - (n + 1)):(df.at[i, 'POS'] + n)].seq
        return(df)


    @staticmethod
    def replacement(string, pos, char):
        '''Replace a substring within a string (e.g. to substitute DNA bases)'''
        string = string[0:pos] + string[pos:(pos + 1)].replace(string[pos], char) + string[(pos + 1):]
        return(string)


    @staticmethod
    def revcomp(seq):
        '''Reverse complement a DNA sequence'''
        comp = str.maketrans('ATCG', 'TAGC')
        return(seq.upper().translate(comp)[::-1])


    @staticmethod
    def reverse_strand(strand):
        '''Reverse the strand of a DNA sequence'''
        reversal = str.maketrans('+-', '-+')
        return(strand.translate(reversal))


    @staticmethod
    def alt_seq(df, n=50):
        '''Create sequences containing the alternative allele'''
        df2 = pd.DataFrame()
        for i in df.index:
            c_A = df[i:(i+1)].copy()
            c_B = df[i:(i+1)].copy()
            c_A.at[i, 'ID'] = c_A.at[i, 'ID'] + '.REF'
            c_B.at[i, 'ID'] = c_B.at[i, 'ID'] + '.ALT'
            if len(df.at[i, 'REF']) > len(df.at[i, 'ALT']): # deletion
                c_B.at[i, 'SEQ'] = df.at[i, 'SEQ'][:n] + df.at[i, 'ALT'] + df.at[i, 'SEQ'][(n + len(df.at[i, 'REF'])):]
            elif len(df.at[i, 'ALT']) > len(df.at[i, 'REF']): # insertion
                c_B.at[i, 'SEQ'] = DataHandler.replacement(df.at[i, 'SEQ'][:-(len(df.at[i, 'ALT']) - 1)], (n), df.at[i, 'ALT'])
            else:
                c_B.at[i, 'SEQ'] = DataHandler.replacement(c_B.at[i, 'SEQ'], n, c_B.at[i, 'ALT'])
            df2 = df2.append([c_A, c_B], ignore_index=True)
        df2 = df2.reset_index(drop=True)
        return(df2)


class BindingHandler:
    def __init__(self):
        pass


    @staticmethod
    def check_variants(variants, pfm_paths, binding=50, pseudocount=0.8, bg={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}, icp=1.0, icb=0.5):
        '''Analyze changes in transcription factor binding based on genomic variants
        binding: <represents the quantile between the minimal and the maximal possible value from the PWM> from searchSeq function in tfbs tools'''
        output = pd.DataFrame()
        for i in pfm_paths:
            tf = i.split('/')[-1].split('_')
            tf_id = tf[0]
            tf_name = tf[1]
            pfm_matrix = pd.read_csv(i, sep='\t', header=None, index_col=0)
            tf_length = pfm_matrix.shape[1]
            pwm_matrix = pfm_matrix.copy()
            for i in ['A', 'C', 'G', 'T']:
                pwm_matrix.loc[i] = pwm_matrix.loc[i] + (bg[i] * pseudocount)
            pwm_matrix = pwm_matrix / pwm_matrix.sum()
            for i in ['A', 'C', 'G', 'T']:
                pwm_matrix.loc[i] = pwm_matrix.loc[i] / bg[i]
            pwm_matrix = np.log2(pwm_matrix)
            tf_max = pwm_matrix.max().sum()
            tf_min = pwm_matrix.min().sum()
            tf_threshold = np.percentile([tf_min, tf_max], binding)
            icm_matrix = pfm_matrix.copy()
            for i in ['A', 'C', 'G', 'T']:
                icm_matrix.loc[i] = icm_matrix.loc[i] + (bg[i] * pseudocount)
            icm_matrix = icm_matrix / icm_matrix.sum()
            icm_matrix_log = np.log2(icm_matrix)
            icm_matrix_sum = (icm_matrix * icm_matrix_log).sum() + 2
            icm_matrix_sum = list(icm_matrix_sum)
            icm_matrix_sum_rev = icm_matrix_sum.copy() # prepare
            icm_matrix_sum_rev.reverse() # reverse
            icm_matrix = icm_matrix * icm_matrix_sum
            icm_cols = icm_matrix.columns # prepare
            icm_matrix_rev = icm_matrix.sort_index(axis=1, ascending=False) # reverse column order to get the correct IC values later on
            icm_matrix_rev.columns = icm_cols # change col names
            for j in variants.index:
                if 'REF' in variants.at[j, 'ID']:
                    if len(variants.at[j, 'REF']) == len(variants.at[j, 'ALT']):
                        indel = False
                    else:
                        indel = True
                    possible_seqs_ref = variants.at[j, 'SEQ'][50 - (tf_length - 1):50 + tf_length] # standard orientation # this means that the first position that is tested is actually the last one in the motif (starting the sequence from the left, the variant of interest will be at the right-most position)
                    if len(variants.at[j, 'REF']) > len(variants.at[j, 'ALT']): # deletion
                        offset = len(variants.at[j, 'REF']) - len(variants.at[j, 'ALT'])
                    else:
                        offset = 0
                    possible_seqs_alt = variants.at[(j + 1), 'SEQ'][50 - (tf_length - 1) - offset:50 + tf_length - offset]
                    scores = pd.DataFrame(columns = ['REF_SCORE', 'ALT_SCORE', 'DIFF', 'RANGE', 'CHANGE', 'REF_THRESHOLD', 'ALT_THRESHOLD', 'CHROM', 'POS', 'REF', 'ALT', 'ID', 'TF_NAME', 'TF_ID', 'TF_LENGTH', 'ORIENTATION', 'ICP', 'ICB_REF', 'ICB_ALT', 'ICP_THRESHOLD', 'ICB_REF_THRESHOLD', 'ICB_ALT_THRESHOLD']) # removed: SILENT, DIFF_PERC, VAR_INDEX, PERC_THRESHOLD
                    for k in range(tf_length): # for every possible position of the variant within the TFBS...
                        current_score_ref = 0
                        current_score_alt = 0
                        if indel:
                            scores.at[k, 'ICP'] = 99 # to pass filtering
                            scores.at[k, 'ICB_REF'] = 99
                            scores.at[k, 'ICB_ALT'] = 99
                        else:
                            scores.at[k, 'ICP'] = icm_matrix_sum_rev[k]
                            scores.at[k, 'ICB_REF'] = icm_matrix_rev.at[variants.at[j, 'REF'], k + 1]
                            scores.at[k, 'ICB_ALT'] = icm_matrix_rev.at[variants.at[j, 'ALT'], k + 1]
                        scores.at[k, 'VAR_POS'] = tf_length - k
                        for l in range(tf_length): # ...calculate the snake sum of the entire motif
                            current_score_ref += pwm_matrix.at[possible_seqs_ref[l + k], l + 1]
                            current_score_alt += pwm_matrix.at[possible_seqs_alt[l + k], l + 1]
                        scores.at[k, 'REF_SCORE'] = current_score_ref
                        scores.at[k, 'ALT_SCORE'] = current_score_alt
                    scores['ICP_THRESHOLD'] = scores['ICP'].apply(lambda x: x >= icp)
                    scores['ICB_REF_THRESHOLD'] = scores['ICB_REF'].apply(lambda x: x >= icb)
                    scores['ICB_ALT_THRESHOLD'] = scores['ICB_ALT'].apply(lambda x: x >= icb)
                    scores['REF_THRESHOLD'] = scores['REF_SCORE'].apply(lambda x: x >= tf_threshold)
                    scores['ALT_THRESHOLD'] = scores['ALT_SCORE'].apply(lambda x: x >= tf_threshold)
                    #scores['VAR_INDEX'] = scores.apply(lambda x: x.index)
                    scores['CHROM'] = variants.at[j, 'CHROM']
                    scores['POS'] = variants.at[j, 'POS']
                    scores['REF'] = variants.at[j, 'REF']
                    scores['ALT'] = variants.at[j, 'ALT']
                    scores['ID'] = variants.at[j, 'ID'].split('.')[0]
                    scores['TF_NAME'] = tf_name
                    scores['TF_ID'] = tf_id
                    scores['TF_LENGTH'] = tf_length
                    scores['ORIENTATION'] = 'standard'
                    scores['RANGE'] = abs(tf_min - tf_max)
                    scores = scores[scores['REF_THRESHOLD'] | scores['ALT_THRESHOLD']]
                    scores = scores[scores['ICP_THRESHOLD']]
                    scores = scores[scores['ICB_REF_THRESHOLD'] | scores['ICB_ALT_THRESHOLD']]
                    output = output.append(scores, ignore_index=True, sort=False)
                    possible_seqs_ref = DataHandler.revcomp(possible_seqs_ref) # revcomp orientation
                    possible_seqs_alt = DataHandler.revcomp(possible_seqs_alt)
                    scores = pd.DataFrame(columns = ['REF_SCORE', 'ALT_SCORE', 'DIFF', 'RANGE', 'CHANGE', 'REF_THRESHOLD', 'ALT_THRESHOLD', 'CHROM', 'POS', 'REF', 'ALT', 'ID', 'TF_NAME', 'TF_ID', 'TF_LENGTH', 'ORIENTATION', 'ICP', 'ICB_REF', 'ICB_ALT', 'ICP_THRESHOLD', 'ICB_REF_THRESHOLD', 'ICB_ALT_THRESHOLD']) # removed: SILENT, DIFF_PERC, VAR_INDEX, PERC_THRESHOLD
                    for k in range(tf_length):
                        current_score_ref = 0
                        current_score_alt = 0
                        if indel:
                            scores.at[k, 'ICP'] = 99
                            scores.at[k, 'ICB_REF'] = 99
                            scores.at[k, 'ICB_ALT'] = 99
                        else:
                            scores.at[k, 'ICP'] = icm_matrix_sum_rev[k]
                            scores.at[k, 'ICB_REF'] = icm_matrix_rev.at[DataHandler.revcomp(variants.at[j, 'REF']), k + 1]
                            scores.at[k, 'ICB_ALT'] = icm_matrix_rev.at[DataHandler.revcomp(variants.at[j, 'ALT']), k + 1]
                        scores.at[k, 'VAR_POS'] = tf_length - k
                        for l in range(tf_length):
                            current_score_ref += pwm_matrix.at[possible_seqs_ref[l + k], l + 1]
                            current_score_alt += pwm_matrix.at[possible_seqs_alt[l + k], l + 1]
                        scores.at[k, 'REF_SCORE'] = current_score_ref
                        scores.at[k, 'ALT_SCORE'] = current_score_alt
                    scores['ICP_THRESHOLD'] = scores['ICP'].apply(lambda x: x >= icp)
                    scores['ICB_REF_THRESHOLD'] = scores['ICB_REF'].apply(lambda x: x >= icb)
                    scores['ICB_ALT_THRESHOLD'] = scores['ICB_ALT'].apply(lambda x: x >= icb)
                    scores['REF_THRESHOLD'] = scores['REF_SCORE'].apply(lambda x: x >= tf_threshold)
                    scores['ALT_THRESHOLD'] = scores['ALT_SCORE'].apply(lambda x: x >= tf_threshold)
                    #scores['VAR_INDEX'] = scores.apply(lambda x: x.index) # not correct?
                    scores['CHROM'] = variants.at[j, 'CHROM']
                    scores['POS'] = variants.at[j, 'POS']
                    scores['REF'] = variants.at[j, 'REF']
                    scores['ALT'] = variants.at[j, 'ALT']
                    scores['ID'] = variants.at[j, 'ID'].split('.')[0]
                    scores['TF_NAME'] = tf_name
                    scores['TF_ID'] = tf_id
                    scores['TF_LENGTH'] = tf_length
                    scores['ORIENTATION'] = 'reverse'
                    scores['RANGE'] = abs(tf_min - tf_max)
                    scores = scores[scores['REF_THRESHOLD'] | scores['ALT_THRESHOLD']]
                    scores = scores[scores['ICP_THRESHOLD']]
                    scores = scores[scores['ICB_REF_THRESHOLD'] | scores['ICB_ALT_THRESHOLD']]
                    output = output.append(scores, ignore_index=True, sort=False)
        return(output)


    @staticmethod
    def filter_variants(df):
        '''Filter output of <check_variants>'''
        if len(df) > 0:
            df['DIFF'] = df.apply(lambda x: x['ALT_SCORE'] - x['REF_SCORE'], axis=1)
            #df['DIFF_PERC'] = df.apply(lambda x: x['DIFF'] / x['REF_SCORE'], axis=1)
            #df['SILENT'] = df['DIFF_PERC'].apply(lambda x: x == 0) # for filtering purposes only
            df['CHANGE'] = df.apply(lambda x: 'GAIN' if x['ALT_SCORE'] > x['REF_SCORE'] else ('LOSS' if x['ALT_SCORE'] < x['REF_SCORE'] else 'SILENT'), axis=1)
            #df['PERC_THRESHOLD'] = df['DIFF_PERC'].apply(lambda x: abs(x) >= min_change)
            #df['PERC_THRESHOLD'] = True # temporarily
            #df['VAR_POS'] = df.apply(lambda x: x['TF_LENGTH'] - x['VAR_INDEX'], axis=1) # correct for minus strand?
            #df = df[df['PERC_THRESHOLD'] | df['SILENT']] # can this lead to an empty df as well?
            df = df.drop(['REF_THRESHOLD', 'ALT_THRESHOLD', 'ICP_THRESHOLD', 'ICB_REF_THRESHOLD', 'ICB_ALT_THRESHOLD'], axis=1) # removed: SILENT, VAR_INDEX
            return(df)
        else: # df is empty, therefore nothing found
            print('No altered transcription factor binding sites found (based on the current parameters).')
            print('Exiting.')
            sys.exit()
