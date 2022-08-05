# -*- coding: utf-8 -*-
"""
Created on Mon May  2 10:38:33 2022

@author: Toño
"""

def block_dict(query_seq, subj_seq, dictionary):
    '''
    Takes two protein sequences (one from the query and the other from the subject), divides them into blocks according to the stops '*' present in both and then saves
    in a list of dictionaries the start, end and score (Matrix Blosum) of each block. Afterwards, the block with better score is selected and the data of the
    score, start, end and sequence from both query and subject is updated in the general dataframe.

    Parameters
    ----------
    query_seq : String
        Protein sequence aligned of the query
    subj_seq : String
        Protein sequence aligned of the subject
    dictionary : Dictionary
        Matrix Blosum Values

    Returns
    -------
    max_score_fragment : 

    '''
    
    list_stops = list()
    
    for index, x in enumerate(query_seq):
        if index == 0 or index == len(query_seq) - 1:
            if x != '*' and subj_seq[index] != '*':
                list_stops.append(index)
                continue
        if x == '*' or subj_seq[index] == '*':
            list_stops.append(index)
        
    block_dict_list = list()
    
    for idx, x in enumerate(list_stops):
        if x == list_stops[-1]:
            continue
        block_dict = dict()
        block_dict['Align_Start'] = x + 1
        block_dict['Align_End'] = list_stops[idx + 1]
        block_dict['Score'] = score(query_seq[x : list_stops[idx + 1]], subj_seq[x : list_stops[idx + 1]], dictionary)
        block_dict_list.append(block_dict)
    
    best_block = sorted(block_dict_list, key = lambda x: x['Score'], reverse=True)
    
    return best_block[0]

def score(query_frag, subj_frag, dictionary):
    '''
    Calculate the score using Matrix Blosum

    Parameters
    ----------
    query_frag : String
        Protein sequence from the query
    subj_frag : String
        Protein sequence from the subject
    dictionary : Dictionary
        Matrix Blosum Values

    Returns
    -------
    Score : Int
        Value according to Matrix Blosum

    '''
    
    score = 0
    for index, x in enumerate(query_frag):
        if x == '-' or subj_frag[index] == '-':
            continue
        score += dictionary[x][subj_frag[index]]
        
    return score
    