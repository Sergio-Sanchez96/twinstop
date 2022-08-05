#!/usr/bin/env -S python -u
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 11:55:40 2022

@author: Toño
"""

import sys
sys.path.insert(0, '/home/mmariotti/bioinfo_dev/pyranges/pyranges')
import pyranges as pr, subprocess, shlex, pandas as pd
from easyterm import command_line_options, check_file_presence, write
from os import path, makedirs, replace, remove
import matplotlib
from matplotlib import pyplot as plt
from easybioinfo import count_coding_changes, count_coding_sites, translate
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from block_selection import block_dict, score
from dotplot_twinstop import dot_plot
from datetime import datetime
import tracemalloc

def run_tblastx(query, db, outfolder, n_cpu, default_file, 
               format_6_file, q_name, s_name, resultados):
    '''
    Runs tblastx in default (optionally) and tabular 6 format.
    (Columns selected: subject accession(sacc), subject title (stitle), 
     subject alignment start (sstart), subject alignment end (send), 
     subject frame (sframe), query accession (qacc), 
     query alignment start (qstart), query alignment end (qend), 
     bitscore (bitscore), evalue (evalue), 
     query aligned protein sequence (qseq), 
     subject aligned protein sequence (sseq)).
    
    Parameters
    ----------
    query : String
        Path where the transcriptome from the query is located.
    db : String
        Path where the transcriptome from the subject is located. It must be 
        previously recognized by Blast with the command 'makeblastdb'.
    outfolder : Easyterm object
        Path to the output folder (opt['o']).
    n_cpu : Int, Easyterm object
        opt['c']
    default_file : Boolean value, easyterm object
        opt['b']
    format_6_file : Boolean value, easyterm object
        opt['f']
    q_name : String
        Name used for the query in the outfile name.
    s_name : String
        Name used for the subject in the outfile name.
    
    Raises
    ------
    Exception
        We raise an exception to stop the program in case returncode returns 
        different than zero, indicating that subprocess.run hasn't run 
        successfully.
        We also print the stdout and stderr to know more about the problem.
    
    Returns
    -------
    default_outfile: Tblastx outfile
        With the results of the tblastx hits in its default version.
        (This file is saved but not returned)
    format_6_outfile : Tblastx outfile
        With the results of the tblastx hits in the tabular format 6.
        
    '''
    
    if default_file == True: # decide about the creation of the default file.
    
        write(f'Running default tblastx')
        
        # name and location of the default file.
        temp_default_outfile = (outfolder + 'temp_default_tblastx_hits.txt')
        default_outfile = (outfolder + 'default_tblastx_hits.txt')
        # string to run tblastx with the query, subject, outfile and number of
        # cpu used.
        default_cmd = ('tblastx -query ' + query + ' -db ' + db + ' -out ' + 
                       temp_default_outfile + ' -num_threads ' + str(n_cpu))
        # shlex.split, splits the string into a shell-like syntax.
        default_cmd_list = shlex.split(default_cmd)
        # subprocess module allows you to spawn new processes, connect to 
        # their input/output/error pies, and obtain their return codes.
        x = subprocess.run(default_cmd_list, capture_output=True)
        
        if x.returncode != 0: 
            print(x.stderr, x.stdout)
            raise Exception()
        # exit status of the process. Typically, an exit status of 0 indicates 
        # that it ran successfully.
        # in case the program fails, we print the standard output (stdout) and
        # the standard error (stderr).
        # an exception is raised to stop the program.
        
        # replace the temp file to be sure than the tblastx is completed
        replace(temp_default_outfile, default_outfile)
        
        write(f'Default tblastx ran successfully, results in {default_outfile}')
    
    # we run the tblastx format 6 table only if it does not exist or if we 
    # force it.
    if path.exists(resultados) == False or format_6_file == True:
        
        write(f'Running format 6 tblastx')
        
        temp_format_6_outfile = outfolder + 'temp_tblastx.tsv'
        # command to run tblastx with the specific columns we want
        format_6_cmd = (
            'tblastx -query ' + query + ' -db ' + db + 
            " -outfmt '6 sacc sstart send sframe qacc qstart" + 
            " qend qframe bitscore evalue sseq qseq' -out " + 
            temp_format_6_outfile + ' -num_threads ' + str(n_cpu))
        format_6_cmd_list = shlex.split(format_6_cmd)
        y = subprocess.run(format_6_cmd_list, capture_output=True)
        
        if y.returncode != 0:
            print(y.stderr, y.stdout)
            raise Exception()
        
        replace(temp_format_6_outfile, resultados)
        
        write(f'Format 6 tblastx ran successfully, results in {resultados}')
    
    return resultados

def conv_pandas_df(format_6_tblastx_alignments):
    '''
    Transforms the tblastx table in two pandas DataFrames, subject and query.
    
    Parameters
    ----------
    format_6_tblastx_alignments : Tblastx tabular format 6 table
        With all the tblastx hits.
    
    Returns
    -------
    table_query : Dataframe
        Dataframe only with the query-related columns.
    table_subj : Dataframe
        Dataframe only with the subject-related columns.
        
    '''
    
    write(f'Converting {format_6_tblastx_alignments} into a pandas DataFrame')
    # creates a DataFrame and we name the columns.
    table_subj = pd.read_table(
        format_6_tblastx_alignments, 
        names=['Chromosome', 'Start', 'End', 'Subj_fr', 'Q_ID', 
               'Q_align_s', 'Q_align_e', 'Q_fr', 'Score', 'Evalue', 
               'Subj_align_prot_seq', 'Q_align_prot_seq'])
    # where 'Strand' is greater than 0 it will be replace by '+',
    # if not, by '-'. Pyranges format.
    table_subj['Strand'] = (
        (table_subj['Subj_fr'] > 0).replace(
            {True:'+', False:'-'}))
    # creates a boolean Series
    indexer = table_subj['Start'] > table_subj['End']
    indexer2 = table_subj['Q_align_s'] > table_subj['Q_align_e']
    # switch the columns of 'Start' and 'End' of both the query and the 
    # subject if they are in a negative frame.
    table_subj.loc[
        indexer, ['Start', 'End']] = table_subj.loc[
            indexer, ['End', 'Start']].values
    # substitutes the value of the column only where indexer is True.
    table_subj.loc[
        indexer2, ['Q_align_s', 'Q_align_e']] = table_subj.loc[
            indexer2, ['Q_align_e', 'Q_align_s']].values
    
    table_subj['ID'] = range(
        1, len(table_subj) + 1) # Identification column (IMPORTANT!!!)
    # BLAST is a 1-base program (does not take the last value), 
    # while Pyranges is 0-base (takes the first and last value).
    table_subj['Start'] = table_subj['Start'] - 1
    table_subj['Q_align_s'] = table_subj['Q_align_s'] - 1
    # divide the dataframe into two dataframes (one for the subject and the 
    # other for the query), ready to be transformed into Pyranges.
    
    return query_subject_dfs(table_subj) # we return both dataframes.

def query_subject_dfs(table_subj):
    '''
    Divides the whole dataframe into two, one for the query and the 
    other for the subject, with the format needed to convert them into 
    Pyranges objects.

    Parameters
    ----------
    table_subj : Dataframe
        Dataframe with all the columns about the blasthits.

    Returns
    -------
    table_query : Dataframe
        Dataframe only with the query-related columns.
    table_subj : Dataframe
        Dataframe only with the subject-related columns.
        
    '''
    
    write(f'Dividing columns into query and subject dataframes')
    
    table_query = table_subj.copy() # creates table_query from the table_subj.
    # drops the query-related columns.
    table_subj.drop(['Q_ID', 'Q_align_s', 'Q_align_e', 'Q_fr',
                     'Q_align_prot_seq'], axis=1, inplace=True)
    # renames the columns to fit into the Pyranges format (Chromosome, Start, 
    # End, Strand (+/-)).
    table_query = table_query.rename(
        columns={'Chromosome': 'Subj_ID', 'Start': 'Subj_align_s', 
                 'End': 'Subj_align_e', 'Q_ID': 'Chromosome', 
                 'Q_align_s': 'Start', 'Q_align_e': 'End'})
    # drops the subject-related columns ('Score' and 'Evalue' are left in 
    # the subj_table).
    table_query.drop(['Strand', 'Subj_ID', 'Subj_align_s', 
                      'Subj_align_e', 'Subj_fr',
                      'Subj_align_prot_seq', 'Score', 'Evalue'], 
                     axis=1, inplace=True)
    table_query['Strand'] = table_query['Q_fr'].copy()
    # Strand column needs to have '+' or '-' only.
    table_query['Strand'] = (
        table_query['Strand'] > 0).replace({True:'+', False:'-'})
    # we return both the query and subject dataframes.
    
    write(f'Returning query and subject dataframes')
    
    return table_query, table_subj

def conv_pyran(tblastx_hits, db_file, query_file):
    '''
    Converts dataframes into pyranges, extract the CDS and protein sequences 
    from both query and subject, and then rejoins them into one (after 
    transforming into dataframes again).
    
    Parameters
    ----------
    tblastx_hits : Dataframe
        Dataframe with the tblastx hits.
    db_file : String
        Path where the transcriptome from the subject is located.
    query_file : String
        Path where the transcriptome from the query is located.
    
    Returns
    -------
    final_table_df : Dataframe
        Dataframe table with the tblastx columns plus the results of the 
        pyranges CDS sequences.
        
    '''
    
    write(f'Converting subject and query dataframes into pyranges')
    
    query_table, subject_table = conv_pandas_df(tblastx_hits)
    
    # control to check that tblastx and pyranges protein sequences are the 
    # same.
    # query_table, subj_table = translate_prot(subject_table, query_table, 
                                             #db_file, query_file, 
                                             #genetic_code = False)
    query_table, subj_table = translate_prot(subject_table, query_table, 
                                             db_file, query_file)
    
    query_table.drop(['Strand'], axis=1, inplace=True)
    # is needed to rename query's columns before joining the two dataframes
    query_table = query_table.rename(columns={'Chromosome': 'Q_ID', 
                                              'Start': 'Q_align_s', 
                                              'End': 'Q_align_e'})
    # joins subject and query DataFrames according to ID columns shared values           
    final_table_df = subject_table.join(query_table.set_index('ID'), on='ID')
    final_table_df = final_table_df.reindex(
        columns=['ID', 'Chromosome', 'Start', 'End', 'Strand', 'Subj_fr', 
                 'Q_ID','Q_align_s', 'Q_align_e', 'Q_fr', 'Score', 
                 'Evalue', 'Subj_CDS', 'Query_CDS', 'Subj_align_prot_seq', 
                 'Q_align_prot_seq'])
    
    return final_table_df # returns the joined dataframe with the new index

def max_score_block(selected_IDs, dictionary_matrix):
    '''
    Selection of the fragment (between two stops) with best score.

    Parameters
    ----------
    selected_IDs : Dataframe
        Blast hits dataframe.
    dictionary_matrix : Dictionary of dictionaries
        Dictionary with the Matrix BLOSUM62 values.

    Returns
    -------
    selected_IDs : Dataframe
        Fragments with the best score from every tblastx-hit.
        
    '''
    
    # lists to update the general dataframe.
    list_subj_start = list()
    list_subj_end = list()
    list_query_start = list()
    list_query_end = list()
    list_query_prot = list()
    list_subj_prot = list()
    list_score = list()
    
    write(f'Taking the fragments with the best score')
    
    # we iter the dafaframe by rows.
    for i, row in selected_IDs.iterrows():
        # selection of the fragment with higher score.
        max_score_frag = block_dict(row['Q_align_prot_seq'], 
                                    row['Subj_align_prot_seq'], 
                                    dictionary_matrix)
        
        list_query_prot.append(row['Q_align_prot_seq'][max_score_frag['Align_Start'] : max_score_frag['Align_End']])
        list_subj_prot.append(row['Subj_align_prot_seq'][max_score_frag['Align_Start'] : max_score_frag['Align_End']])
        list_score.append(max_score_frag['Score'])
        # for the positive frames.
        if row['Subj_fr'] > 0:
            # if the fragment is at the beginning of the sequence
            if max_score_frag['Align_Start'] == 1:
                row['End'] = (
                    row['End'] - 
                    (len(row['Subj_align_prot_seq']) -
                     max_score_frag['Align_End']) * 3)
            # if the fragment is at the end of the sequence
            elif max_score_frag['Align_End'] == len(row['Subj_align_prot_seq']):
                row['Start'] = (
                    row['Start'] + 
                    (len(row['Subj_align_prot_seq']) - 
                     max_score_frag['Align_End']) * 3)
            # if the fragment is inbetween the sequence
            else:
                row['Start'] = (
                    row['Start'] + 
                    max_score_frag['Align_Start'] * 3)
                row['End'] = (
                    row['End'] - 
                    (len(row['Subj_align_prot_seq']) - 
                     max_score_frag['Align_End']) * 3)
        else: # negative frames
            if max_score_frag['Align_Start'] == 1:
                row['Start'] = (
                    row['Start'] + 
                    (len(row['Subj_align_prot_seq']) - 
                     max_score_frag['Align_End']) * 3)
            elif max_score_frag['Align_End'] == len(row['Subj_align_prot_seq']):
                row['End'] = (
                    row['End'] - 
                    (len(row['Subj_align_prot_seq']) - 
                     max_score_frag['Align_End']) * 3)
            else:
                row['Start'] = (
                    row['Start'] + 
                    (len(row['Subj_align_prot_seq']) - 
                     max_score_frag['Align_End']) * 3)
                row['End'] = (
                    row['End'] - 
                    max_score_frag['Align_Start'] * 3)
                
        list_subj_start.append(row['Start'])
        list_subj_end.append(row['End'])
        
        if row['Q_fr'] > 0:
            if max_score_frag['Align_Start'] == 1:
                row['Q_align_e'] = (
                    row['Q_align_e'] - 
                    (len(row['Q_align_prot_seq']) - 
                     max_score_frag['Align_End']) * 3)
            elif max_score_frag['Align_End'] == len(row['Q_align_prot_seq']):
                row['Q_align_s'] = (
                    row['Q_align_s'] + 
                    (len(row['Q_align_prot_seq']) - 
                     max_score_frag['Align_End']) * 3)
            else:
                row['Q_align_s'] = (
                    row['Q_align_s'] + 
                    max_score_frag['Align_Start'] * 3)
                row['Q_align_e'] = (
                    row['Q_align_e'] - 
                    (len(row['Q_align_prot_seq']) - 
                     max_score_frag['Align_End']) * 3)
        else:
            if max_score_frag['Align_Start'] == 1:
                row['Q_align_s'] = (
                    row['Q_align_s'] + 
                    (len(row['Q_align_prot_seq']) - 
                     max_score_frag['Align_End']) * 3)
            elif max_score_frag['Align_End'] == len(row['Q_align_prot_seq']):
                row['Q_align_e'] = (
                    row['Q_align_e'] - 
                    (len(row['Q_align_prot_seq']) - 
                     max_score_frag['Align_End']) * 3)
            else:
                row['Q_align_s'] = (
                    row['Q_align_s'] + 
                    (len(row['Q_align_prot_seq']) - 
                     max_score_frag['Align_End']) * 3)
                row['Q_align_e'] = (
                    row['Q_align_e'] - 
                    max_score_frag['Align_Start'] * 3)
                
        list_query_start.append(row['Q_align_s'])
        list_query_end.append(row['Q_align_e'])
                
    # updates the dataframe columns
    selected_IDs['Start'] = list_subj_start
    selected_IDs['End'] = list_subj_end
    selected_IDs['Q_align_s'] = list_query_start
    selected_IDs['Q_align_e'] = list_query_end
    selected_IDs['Q_align_prot_seq'] = list_query_prot
    selected_IDs['Subj_align_prot_seq'] = list_subj_prot
    selected_IDs['Score'] = list_score
    
    write(f'Updated dataframe with the best fragments')
    
    return selected_IDs

def dictionary_seleno():
    '''
    Creates a dictionary of dictionaries with the values of the BLOSUM Matrix.
    
    Returns
    -------
    dictionary_sel : Dictionary
        Dictionary of dictionaries with the values of the BLOSUM Matrix.
        
    '''
    
    write(f'Writing the dictionary BLOSUM62')
    
    dictionary_sel = dict()
    with open('/home/ssanchez/seleno_prediction/enlaces_fungi/Matrix_BLOSUM62sel.txt', 'r') as fr:
        for index, row in enumerate(fr):
            spt = row.split(' ')
            # filter None values
            new_spt = list(filter(None, spt))
            if index == 0:
                # delete empthy spaces
                header = [x.strip() for x in new_spt] 
            spt = [x.strip() for x in new_spt]
            new_spt = list(filter(None, spt))
            if index > 0:
                # converts the values (string) into integers
                ints = [int(x) for x in new_spt[1:]]
                n = 0
                for idx, column in enumerate(header):
                    if n == 0:
                        # create different dictionaries inside one
                        dictionary_sel[header[index-1]] = {}
                        n += 1
                    dictionary_sel[header[index-1]][column] = ints[idx]
                    
    return dictionary_sel

def overlapping(real_final_table):
    '''
    First filter of the script.
    
    Based on selecting only the best score among the overlapping hits.
    
    Parameters
    ----------
    real_final_table : Dataframe
        With all the columns.
    
    Returns
    -------
    selected_IDs : Dataframe
        Dataframe table with the best hit of the overlapping sequences.
        
    '''
    
    write(f'Overlapping')
    # converts into PyRanges
    final_table_pr = pr.PyRanges(real_final_table)
    # creates a 'Cluster' column identifying the overlapping sequences
    final_table_pr = final_table_pr.cluster(strand = False, slack = 0)
    # converts into Dataframe
    final_table_df = final_table_pr.as_df()
    # discards the evalues greater than 0.05
    final_table_df = final_table_df[final_table_df['Evalue'] < 0.05]
    # creates a series with as many rows as clusters in the dataframe
    cluster_ID = final_table_df['Cluster'].unique()
    # creates an empthy dataFrame
    selected_IDs = pd.DataFrame()
    
    for i in cluster_ID:
        # keeps only those rows where 'Cluster' column is equal to 'i'
        cluster = final_table_df[final_table_df['Cluster'] == i]
        # keeps the row with the best score
        cluster = cluster[cluster['Score'] == cluster['Score'].max()]
        # appends the selected row in selected_IDs dataframe, to ignore_index 
        # is needed.
        selected_IDs = selected_IDs.append(cluster, ignore_index=True)
    
    # returns the new Dataframe with only the rows with the best scores 
    # among the overlapping hits.
    return selected_IDs

def conv_gff(clusters, query_gff, subject_gff):
    '''
    Converts dataframes into gffs.
    
    Parameters
    ----------
    clusters : Dataframe
        General Dataframe with all the columns.
    query_gff : String
        Path where the gff of the query will be saved.
    subject_gff : String
        Path where the gff of the subject will be saved.
    
    Returns
    -------
    None
    
    '''
    
    write(f'GFF format')
    
    table_query, table_subj = query_subject_dfs(clusters)
    # we have to change 'ID' column's name to fit with the gff format.
    table_query = table_query.rename(columns={'ID':'Attribute'})
    # we also change 'Score' column's name to include it in the 'Attribute' 
    # column.
    # 'extend_orfs.py' does not take the score values for the 'Score' column
    table_subj = table_subj.rename(columns={'ID':'Attribute', 
                                            'Score':'Value'}) 
    # creates a new column called 'Feature' with the same value ('CDS') 
    # in all rows.
    table_query['Feature'] = 'CDS'
    table_subj['Feature'] = 'CDS'
    # we only take the columns we need to fulfill the gff format
    table_query_reduced = table_query.loc[:, ['Chromosome', 'Start', 'End', 
                                              'Strand', 'Attribute', 
                                              'Feature', 'Q_fr']]
    # 'Attribute', 'Value' and 'Evalue' will be put together in the 
    # 'Attribute' column.
    table_subj_reduced = table_subj.loc[:, ['Chromosome', 'Start', 'End', 
                                            'Strand', 'Attribute', 'Feature', 
                                            'Value', 'Evalue', 'Subj_fr']]
    
    py_query = pr.PyRanges(table_query_reduced) # converts into PyRanges.
    py_subj = pr.PyRanges(table_subj_reduced)
    
    py_query.to_gff3(path = query_gff) # converts pyranges object into gff.
    py_subj.to_gff3(path = subject_gff)

def extend_orfs(clusters, subject_gff, query_gff, subject, query, 
                out_subj_gff, out_query_gff, output):
    '''
    Runs Marco's script 'extend_orfs.py' which extends the CDS sequences from 
    the query and the subject both downstream and upstream, until a stop is 
    found or the transcript finishes (downstream) and a initial codon is found 
    (Methionine) or the transcript finishes (upstream).
    
    Parameters
    ----------
    clusters : Dataframe
        Dataframe with the best scored tblastx hits among the overlapping ones.
    subject_gff : String
        Path where the gff file related to the subject is saved.
    query_gff : String
        Path where the gff file related to the query is saved.
    subject : String
        Path where the transcriptome of the subject is located.
    query : String
        Path where the transcriptome of the query is located.
    out_subj_gff : String
        Path where the gff file related to the subject will be saved.
    out_query_gff : String
        Path where the gff file related to the query will be saved.
    output : Easyterm object
        opt['o']
        
    Raises
    ------
    Exception
        We raise an exception to stop the program in case returncode returns 
        different than zero, indicating that subprocess.run hasn't run 
        successfully.
        We also print the stdout and stderr to know more about the problem.
        
    Returns
    -------
    joined_df : Dataframe
        Dataframe with the alignments after extending orfs.
        
    '''
    
    write(f'Extending orfs')
    
    conv_gff(clusters, query_gff, subject_gff)
    
    # run 'extend_orfs', one per subject and one per query.
    # the extension stops only when it encounters an 'TAG' or 'TAA' stop.
    subj = ('/home/mmariotti/scripts/extend_orfs.py -i ' + subject_gff + 
            ' -t ' + subject + ' -o ' + out_subj_gff + " -stops 'TAG,TAA'")
    subj_list = shlex.split(subj)
    x = subprocess.run(subj_list, capture_output=True)
    
    if x.returncode != 0:
        print(x.stderr, x.stdout)
        raise Exception()
    
    q = ('/home/mmariotti/scripts/extend_orfs.py -i ' + query_gff + ' -t ' + 
         query + ' -o ' + out_query_gff + " -stops 'TAG,TAA'")
    q_list = shlex.split(q)
    y = subprocess.run(q_list, capture_output=True)
    if y.returncode != 0:
        print(y.stderr, y.stdout)
        raise Exception()
    
    remove(subject_gff) # removes the temporal files.
    remove(query_gff)
    
    # name the columns of the gff resulting file converting it into a 
    # dataframe.
    subj_df_gff = pd.read_csv(out_subj_gff, 
                              sep='\t', names=['Chromosome', 'Source', 
                                               'Feature', 'Start', 'End', 
                                               'Score', 'Strand', 'Frame', 
                                               'Attribute'])
    # All the columns that do not have one of this names will be put together
    # in the 'Attribute' column.
    query_df_gff = pd.read_csv(out_query_gff, 
                               sep='\t', names=['Chromosome', 'Source', 
                                                'Feature', 'Start', 'End', 
                                                'Score', 'Strand', 'Frame', 
                                                'Attribute'])
    remove(out_subj_gff)
    remove(out_query_gff)
    
    query_df_gff['Start'] = query_df_gff['Start'] - 1
    subj_df_gff['Start'] = subj_df_gff['Start'] - 1
    
    return pandas_2(subj_df_gff, query_df_gff, subject, query, output)
    
def translate_prot(subj_df, query_df, subject, query, genetic_code = True):
    '''
    Function to translate the nucleotidic sequences got by PyRanges into 
    protein using translate(), from Marco's easyterm module.
    
    Parameters
    ----------
    subj_df : Dataframe
        Dataframe with the subject-related columns.
    query_df : Dataframe
        Dataframe with the query-related columns.
    subject : String
        Path to the file of the subject.
    query : String
        Path to the file of the query.
    genetic_code : Boolean
        Controls if the protein sequences have 'U's or not.
        
    Returns
    -------
    query_df : Dataframe
        Dataframe  with the query-related columns after translation.
    subj_df : Dataframe
        Dataframe with the subject-related columns after translation.
        
    '''
    
    write(f'Translating into protein')
    
    # convert into gff or extend_orfs changes the 'Start' values to Blast 
    # format, so it is neccessary to subtract one again before transforming 
    # into PyRanges.
    subj_pr = pr.PyRanges(subj_df)
    query_pr = pr.PyRanges(query_df)
    # get the CDS sequences (nucleotides)
    seq_query = pr.get_fasta(query_pr, query)
    seq_subj = pr.get_fasta(subj_pr, subject)
    
    if genetic_code == True:
        query_pr.Query_CDS = seq_query
        # translates the CDS sequences into protein (conserves the 'U's)
        prot_query = [translate(s, genetic_code='1+U') for s in seq_query]
        query_pr.Q_align_prot_seq = prot_query
        
        subj_pr.Subj_CDS = seq_subj
        prot_subj = [translate(s, genetic_code='1+U') for s in seq_subj]
        subj_pr.Subj_align_prot_seq = prot_subj
        
        write(f'Protein sequences with Selenocysteine (U)')
        
    else:
        write(f'Checking protein sequences')
        # translates the CDS sequences into protein (no 'U's)
        prot_query = pd.Series([translate(s) for s in seq_query])
        prot_query.index = query_pr.Q_align_prot_seq.index.copy()
        if prot_query.equals(query_pr.Q_align_prot_seq) == True:
            write(f'Query protein sequences checked')
        else:
            # '~' symbol inverse the characters of a series, True is 
            # transformed in False and the other way around.
            # index = those rows where the 'Q_align_prot_seq' is different
            # (False) to the prot_query.
            index = ~(prot_query.eq(query_pr.Q_align_prot_seq))
            print(seq_query[index])
            print(prot_query[index])
            print(query_pr.Q_align_prot_seq[index])
            raise Exception()
        
        prot_subj = pd.Series([translate(s) for s in seq_subj])
        prot_subj.index = subj_pr.Subj_align_prot_seq.index.copy()
        if prot_subj.equals(subj_pr.Subj_align_prot_seq) == True:
            write(f'Subject protein sequences checked')
        else:
            index = ~(prot_subj.eq(subj_pr.Subj_align_prot_seq))
            print(seq_subj[index])
            print(prot_subj[index])
            print(subj_pr.Subj_align_prot_seq[index])
            raise Exception()
        write(f'Tblastx and Pyranges give the same protein sequences')
        
    query_df = query_pr.as_df()
    subj_df = subj_pr.as_df()
    
    return query_df, subj_df
    
def pandas_2(subj_df_gff, query_df_gff, subject, query, output):
    '''
    This function joins the subject and query gff-format-dataframes.
    
    Parameters
    ----------
    
    subj_df_gff : Dataframe
        Dataframe with the gff columns related to the subject.
    query_df_gff : Dataframe
        Dataframe with the gff columns related to the query.
    subject : String
        Path to the file of the subject.
    query : String
        Path to the file of the query.
    output : Easyterm object
        Path to the output folder.
        
    Returns
    -------
    joined_df : Dataframe
        Dataframe with all the columns of the tblastx hits
        
    '''
    
    query_df, subj_df = translate_prot(subj_df_gff, query_df_gff, 
                                       subject, query)
    
    # rename the PyRanges-format to join both DataFrames
    query_df = query_df.rename(columns={'Chromosome':'Q_ID', 
                                        'Start':'Q_align_s', 
                                        'End':'Q_align_e', 
                                        'Strand':'Q_Strand'})
    # 'extend_orfs.py' is made to ignore these columns, so they are empthy
    query_df.drop(['Source', 'Feature', 'Score', 'Frame'], 
                  axis=1, inplace=True)
    subj_df.drop(['Source', 'Feature', 'Score', 'Frame'], 
                 axis=1, inplace=True)
    
    list_attribute_query = list()
    list_fr_query = list()
    list_attribute_subj = list()
    list_value_subj = list()
    list_fr_subj = list()
    list_evalue_subj = list()
    # now we need to separate the different columns in 'Attribute'
    for i in query_df.index:
        attribute_fr_query = query_df.loc[i, 'Attribute'].split(';')
        for x in attribute_fr_query:
            if x.startswith('Attribute'):
                # strip function removes characters or white spaces from the 
                # beginning or end of a string.
                attribute_query = x.split('=')[1].strip()
                list_attribute_query.append(attribute_query)
            elif x.startswith('Q_fr'):
                fr_query = x.split('=')[1].strip()
                list_fr_query.append(fr_query)
        attribute_value_evalue_fr_subj = subj_df.loc[i, 'Attribute'].split(';')
        for x in attribute_value_evalue_fr_subj:
            if x.startswith('Attribute'):
                attribute_subj = x.split('=')[1].strip()
                list_attribute_subj.append(attribute_subj)
            elif x.startswith('Value'):
                value_subj = x.split('=')[1].strip()
                list_value_subj.append(value_subj)
            elif x.startswith('Subj_fr'):
                fr_subj = x.split('=')[1].strip()
                list_fr_subj.append(fr_subj)
            else:
                evalue_subj = x.split('=')[1].strip()
                list_evalue_subj.append(evalue_subj)
                
    query_df['Attribute'] = list_attribute_query
    query_df['Q_fr'] = list_fr_query
    subj_df['Attribute'] = list_attribute_subj
    subj_df['Score'] = list_value_subj
    subj_df['Subj_fr'] = list_fr_subj
    subj_df['Evalue'] = list_evalue_subj
    
    # merge both dataframes according to 'Attribute' column
    joined_df = subj_df.merge(query_df.set_index(['Attribute']), 
                              on='Attribute')
    return joined_df

def pairwise(joined_df, matrix):
    '''
    This function runs pairwise alignment tool to introduce gaps in the 
    protein sequences.

    Parameters
    ----------
    joined_df : Dataframe
        Dataframe with all the remaining tblastx hits.
    matrix : Dictionary of dictionaries
        Dictionary with the BLOSUM matrix values.

    Returns
    -------
    joined_df : Dataframe
        Dataframe with all the remaining tblastx hits plus gaps in the protein
        sequences.

    '''
    
    write(f'Pairwise alignment')
    # matrix_imp is a dictionary of tuples, with the BLOSUM matrix values
    matrix_imp = matlist.blosum62
    list_pairwise_query = list()
    list_pairwise_subject = list()
    
    for x in matrix:
        for y in matrix[x]:
            if (x, y) in matrix_imp or (y, x) in matrix_imp:
                continue
            else:
                # complete the BLOSUM dictionary with the modified values (U)
                matrix_imp[(x, y)] = matrix[x][y]
                matrix_imp[(y, x)] = matrix[y][x]
                
    for i, row in joined_df.iterrows():
        # -7 is the cost to open a gap, -1 is the cost to extend it
        alignment = pairwise2.align.globalds(row['Q_align_prot_seq'], 
                                             row['Subj_align_prot_seq'], 
                                             matrix_imp, -7, -1,
                                             one_alignment_only=True)
        for x in alignment:
            for idx, y in enumerate(x):
                if idx == 0:
                    list_pairwise_query.append(y)
                elif idx == 1:
                    list_pairwise_subject.append(y)
                else:
                    continue
                
    joined_df['Q_align_prot_seq'] = list_pairwise_query
    joined_df['Subj_align_prot_seq'] = list_pairwise_subject
    
    return joined_df

def UGA(query_prot_seq, subj_prot_seq, u_subj, u_query):
    '''
    Function made to find the 'U' responsible for the read-through.
    
    Parameters
    ----------
    query_prot_seq : String
        Protein sequence corresponding to the query.
    subj_prot_seq : String
        Protein sequence corresponding to the subject.
    u_subj : Int
        Index of the first 'U' in the subject protein sequence.
    u_query : Int
        Index of the first 'U' in the query protein sequence.
        
    Returns
    -------
    u_subj : Int
        Index of the good 'U' in the subject protein sequence.
    u_query : Int
        Index of the good 'U' in the query protein sequence.
        
    '''
    
    list_align_ugas = list()
    center_alignment_subj = len(subj_prot_seq)/2
    center_alignment_q = len(query_prot_seq)/2
    
    for idx, x in enumerate(subj_prot_seq):
        if x == 'U' and query_prot_seq[idx] == 'U':
            list_align_ugas.append(idx)
    # maybe the 'U's are not aligned.
    if len(list_align_ugas) != 0:
        # we choose the aligned 'U' closer to the center of the longer sequence
        if center_alignment_subj >= center_alignment_q:
            # ascending list by default
            # lambda allows to create a function in one line, input : output
            sorted_list = sorted(
                list_align_ugas, 
                key = lambda x: abs(center_alignment_subj - x))
        else:
            sorted_list = sorted(
                list_align_ugas, 
                key = lambda x: abs(center_alignment_q - x))
            
        u_subj = sorted_list[0]
        u_query = sorted_list[0]
    
    return u_subj, u_query

def list_UGAs(row, u_subj, u_query):
    '''
    This function creates 4 lists to separate the U's in upstream or downstream
    according to their position regarding the good selenocysteine (u_subj and 
                                                                   u_query).

    Parameters
    ----------
    row : Row of a Dataframe
        Tblastx hit with all the columns.
    u_subj : Int
        Index of the good 'U' in the subject protein sequence.
    u_query : Int
        Index of the good 'U' in the query protein sequence.

    Returns
    -------
    list_up_query : List
        List of U's indexes upstream the selenocysteine (u_query).
    list_down_query : List
        List of U's indexes downstream the selenocysteine (u_query).
    list_up_subj : List
        List of U's indexes upstream the selenocysteine (u_subj).
    list_down_subj : List
        List of U's indexes downstream the selenocysteine (u_subj).

    '''
    
    list_up_query = list()
    list_down_query = list()
    list_up_subj = list()
    list_down_subj = list()
    
    for idx, x in enumerate(row['Q_align_prot_seq']):
        if x == 'U':
            if idx < u_query:
                list_up_query.append(idx)
            elif idx > u_query:
                list_down_query.append(idx)
        if row['Subj_align_prot_seq'][idx] == 'U':
            if idx < u_subj:
                list_up_subj.append(idx)
            elif idx > u_subj:
                list_down_subj.append(idx)
    
    return list_up_query, list_down_query, list_up_subj, list_down_subj

def UGA_alignments(best_hits, dictionary_matrix, conservation_up, 
                   conservation_down):
    '''
    Second filter of the script.
    
    During this filter, we will get only those hits with aligned selenocysteines
    (U) in query and target protein sequences.
    
    Parameters
    ----------
    best_hits : Dataframe
        Dataframe with all the remaining hits.
    dictionary_matrix : Dictionary of dictionaries
        Dictionary with the BLOSUM matrix values.
    conservation_up : Int, easyterm object
        opt['cons_up']
    conservation_down : Int, easyterm object
        opt['cons_down']
    
    Returns
    -------
    selenoproteins : Dataframe
        Dataframe with the selenoproteins candidates.
        
    '''
    
    write(f'Alignment of UGAs')
    
    list_start_query = list()
    list_stop_query = list()
    list_start_subj = list()
    list_stop_subj = list()
    list_prot_query = list()
    list_prot_subj = list()
    list_cds_query = list()
    list_cds_subj = list()
    list_score = list()
    list_conservation_up = list()
    list_conservation_down = list()
    list_density_score = list()
    
    selenoproteins = pd.DataFrame()
    
    for i, row in best_hits.iterrows():
        # counts the number of U's in both protein sequences
        n_stops_subj = row['Subj_align_prot_seq'].count('U')
        n_stops_query = row['Q_align_prot_seq'].count('U')
        
        if n_stops_subj >= 1 and n_stops_query >= 1:
            # tells the index of the first U in both protein sequences
            u_subj = row['Subj_align_prot_seq'].index('U')
            u_query = row['Q_align_prot_seq'].index('U')
            
            if n_stops_subj >= 2 or n_stops_query >= 2:
                # if there are no aligned U's this function will return the 
                # same u_subj, u_query obtained before with the index('U').
                u_subj, u_query = UGA(row['Q_align_prot_seq'], 
                                      row['Subj_align_prot_seq'],
                                      u_subj, u_query)
            if u_subj == u_query:
                if n_stops_subj >= 2 or n_stops_query >= 2:
                    (list_up_query, list_down_query, 
                     list_up_subj, list_down_subj) = list_UGAs(row, u_subj, 
                                                               u_query)
                    # updates protein and cds sequences according to whether 
                    # there is U's upstream, downstream or both.
                    if len(list_down_query) != 0:
                        row['Q_align_prot_seq'] = (
                            row['Q_align_prot_seq'][:list_down_query[0] + 1] + 
                            '-' * (len(row['Q_align_prot_seq']) - list_down_query[0] - 1))
                        row['Query_CDS'] = row['Query_CDS'][:(list_down_query[0] + 1 -
                                                              row['Q_align_prot_seq'][:list_down_query[0] + 1].count('-')) * 3]
                    if len(list_up_query) != 0:
                        gaps_up = row['Q_align_prot_seq'][:list_up_query[-1] + 1].count('-')
                        row['Q_align_prot_seq'] = (
                            '-' * (list_up_query[-1] + 1) + 
                            row['Q_align_prot_seq'][list_up_query[-1] + 1:])
                        row['Query_CDS'] = row['Query_CDS'][(list_up_query[-1] + 1 - 
                                                             gaps_up) * 3:]

                    if len(list_down_subj) != 0:
                        row['Subj_align_prot_seq'] = (
                            row['Subj_align_prot_seq'][:list_down_subj[0] + 1] + 
                            '-' * (len(row['Subj_align_prot_seq']) - list_down_subj[0] - 1))
                        row['Subj_CDS'] = row['Subj_CDS'][:(list_down_subj[0] + 1 -
                                                           row['Subj_align_prot_seq'][:list_down_subj[0] + 1].count('-')) * 3]
                    if len(list_up_subj) != 0:
                        gaps_up = row['Subj_align_prot_seq'][:list_up_subj[-1] + 1].count('-')
                        row['Subj_align_prot_seq'] = (
                            '-' * (list_up_subj[-1] + 1) +
                            row['Subj_align_prot_seq'][list_up_subj[-1] + 1:])
                        row['Subj_CDS'] = row['Subj_CDS'][(list_up_subj[-1] + 1 -
                                                           gaps_up) * 3:]  
                    # updates the 'Start'/'Q_align_s' or the 'End'/'Q_align_e'
                    # according to the 'Strand'/'Q_strand', respectively.
                    if len(list_up_subj) != 0:
                        if row['Strand'] == '+':
                            row['Start'] = (
                                row['Start'] + 
                                len(row['Subj_align_prot_seq'][:list_up_subj[-1] + 1]) - 
                                row['Subj_align_prot_seq'][:list_up_subj[-1] + 1].count('-'))
                        else:
                            row['End'] = (
                                row['End'] + 
                                len(row['Subj_align_prot_seq'][:list_up_subj[-1] + 1]) - 
                                row['Subj_align_prot_seq'][:list_up_subj[-1] + 1].count('-'))
                    if len(list_down_subj) != 0:
                        if row['Strand'] == '+':
                            row['End'] = (
                                row['End'] - 
                                len(row['Subj_align_prot_seq'][list_down_subj[0] + 1:].replace('-', '')))
                        else:
                            row['Start'] = (
                                row['Start'] - 
                                len(row['Subj_align_prot_seq'][list_down_subj[0] + 1:].replace('-', '')))
                    
                    if len(list_up_query) != 0:
                        if row['Q_Strand'] == '+':
                            row['Q_align_s'] = (
                                row['Q_align_s'] + 
                                len(row['Q_align_prot_seq'][:list_up_query[-1] + 1]) - 
                                row['Q_align_prot_seq'][:list_up_query[-1] + 1].count('-'))
                        else:
                            row['Q_align_e'] = (
                                row['Q_align_e'] + 
                                len(row['Q_align_prot_seq'][:list_up_query[-1] + 1]) - 
                                row['Q_align_prot_seq'][:list_up_query[-1] + 1].count('-'))
                    if len(list_down_query) != 0:
                        if row['Q_Strand'] == '+':
                            row['Q_align_e'] = (
                                row['Q_align_e'] - 
                                len(row['Q_align_prot_seq'][list_down_query[0] + 1:].replace('-', '')))
                        else:
                            row['Q_align_s'] = (
                                row['Q_align_s'] - 
                                len(row['Q_align_prot_seq'][list_down_query[0] + 1:].replace('-', '')))
                    
                    prot_seq_query = ''
                    prot_seq_subj = ''
                    # deletes the non-sense gaps ('-') in both sequences
                    for idx, x in enumerate(row['Subj_align_prot_seq']):
                        if x == '-' and row['Q_align_prot_seq'][idx] == '-':
                            continue
                        else:
                            prot_seq_query += row['Q_align_prot_seq'][idx]
                            prot_seq_subj += x
                    
                    row['Q_align_prot_seq'] = prot_seq_query
                    row['Subj_align_prot_seq'] = prot_seq_subj
                
                n_gaps_subj = row['Subj_align_prot_seq'][:u_subj].count('-')
                n_gaps_query = row['Q_align_prot_seq'][:u_query].count('-')
                # when using cds sequences we need to subtrack the number of
                # gaps and multiply by 3 (1 aa = 3 nucleotides).
                index_3t_nucl_subj = (u_subj - n_gaps_subj) * 3
                index_3t_nucl_query = (u_query - n_gaps_query) * 3
                # filters only the real selenoproteins ('U' = 'TGA')
                if row['Subj_CDS'][index_3t_nucl_subj:index_3t_nucl_subj + 3] == 'TGA' and (
                        row['Query_CDS'][index_3t_nucl_query:index_3t_nucl_query + 3] == 'TGA'):
                    # measures the conservation values of the alignment, according
                    # to the BLOSUM matrix values.
                    conservation_before_tga, conservation_after_tga = (
                        conservation(row, dictionary_matrix, u_query))
                    # filter those hits with values greater than pre-selected
                    # values (conservation_up and conservation_down).
                    if conservation_before_tga >= conservation_up and (
                            conservation_after_tga >= conservation_down):
                        
                        list_conservation_up.append(conservation_before_tga)
                        list_conservation_down.append(conservation_after_tga)
                        list_start_query.append(row['Q_align_s'])
                        list_stop_query.append(row['Q_align_e'])
                        list_start_subj.append(row['Start'])
                        list_stop_subj.append(row['End'])
                        list_prot_query.append(row['Q_align_prot_seq'])
                        list_prot_subj.append(row['Subj_align_prot_seq'])
                        list_cds_query.append(row['Query_CDS'])
                        list_cds_subj.append(row['Subj_CDS'])
                        row['Score'] = score(row['Q_align_prot_seq'], 
                                             row['Subj_align_prot_seq'], 
                                             dictionary_matrix)
                        list_score.append(row['Score'])
                        density_score = round(
                            row['Score']/len(row['Q_align_prot_seq']), 4)
                        list_density_score.append(density_score)
                        # putting the iloc in this way, we print the entire row
                        selenoproteins = (
                            selenoproteins.append(best_hits.iloc[[i]], 
                                                  ignore_index=True))
    
    selenoproteins['Conservation_up'] = list_conservation_up
    selenoproteins['Conservation_down'] = list_conservation_down
    selenoproteins['Q_align_s'] = list_start_query
    selenoproteins['Q_align_e'] = list_stop_query
    selenoproteins['Start'] = list_start_subj
    selenoproteins['End'] = list_stop_subj
    selenoproteins['Q_align_prot_seq'] = list_prot_query
    selenoproteins['Subj_align_prot_seq'] = list_prot_subj
    selenoproteins['Query_CDS'] = list_cds_query
    selenoproteins['Subj_CDS'] = list_cds_subj
    selenoproteins['Score'] = list_score
    selenoproteins['Density Score'] = list_density_score
    
    selenoproteins.sort_values(by='Density Score', inplace=True, 
                               ignore_index=True, ascending=False)
    
    selenoproteins.drop_duplicates(subset='Query_CDS',
                                   inplace=True, ignore_index=True)
    return selenoproteins
    
def conservation(row, dictionary_matrix, u_query):
    '''
    This function measures the score value upstream and downstream of the 
    alignment of query and target.
    
    Parameters
    ----------
    row : Row of a Dataframe
        Tblastx hit with all the columns.
    dictionary_matrix : Dictionary of dictionaries
        Dictionary with the BLOSUM matrix values.
    u_query : Int
        Index of the good 'U' in the query protein sequence.
        
    Returns
    -------
    conservation_before_tga : Int
        Conservation score upstream of the alignment.
    conservation_after_tga : Int
        Conservation score downstream of the alignment.
        
    '''
    
    conservation_before_tga = 0
    conservation_after_tga = 0
    
    for i, x in enumerate(row['Q_align_prot_seq']):
        # gaps and stops are not having into account to measure score
        if x == '*' or x == '-' or row['Subj_align_prot_seq'][i] == (
                '-' or row['Subj_align_prot_seq'][i] == '*'):
            continue
        # score upstream
        if i < u_query:
            conservation_before_tga += dictionary_matrix[x][row['Subj_align_prot_seq'][i]]
        # score downstream
        elif i > u_query:
            conservation_after_tga += dictionary_matrix[x][row['Subj_align_prot_seq'][i]]
         
    return conservation_before_tga, conservation_after_tga

def make_aligned_cds(subj_cds, query_cds, subj_aligned_pep, query_aligned_pep):
    '''
    Function to introduce gaps into the nucleotidic sequences of both query
    and subject.
    
    Parameters
    ----------
    subj_cds : String
        Nucleotidic sequence from the subject.
    query_cds : String
        Nucleotidic sequence from the query.
    subj_aligned_pep : String
        Protein sequence from the subject.
    query_aligned_pep : String
        Protein sequence from the query.
        
    Returns
    -------
    subj_aligned_cds : String
        Nucleotidic sequence from the subject, with gaps.
    query_aligned_cds : String
        Nucleotidic sequence from the query, with gaps.
        
    '''
    
    subj_aligned_cds = subj_cds
    query_aligned_cds = query_cds
    
    for idx, x in enumerate(subj_aligned_pep):
        # for each gap on protein sequence we have to add 3 gaps on the cds 
        # sequence.
        if x == '-':
            subj_aligned_cds = (
                subj_aligned_cds[:idx * 3] + '-' * 3 + 
                subj_aligned_cds[idx * 3:])
    for idx, x in enumerate(query_aligned_pep):
        if x == '-':
            query_aligned_cds = (
                query_aligned_cds[:idx * 3] + '-' * 3 + 
                query_aligned_cds[idx * 3:])
    
    return subj_aligned_cds, query_aligned_cds

def dN_dS(row):
    '''
    This function measures the dN/dS ratio and changes upstream and 
    downstream of the U alignments.
    
    Parameters
    ----------
    row : Row of a Dataframe
        Tblastx hit with all the columns.
        
    Returns
    -------
    u_dN_dS : String/Float
        Non-synonimous/synonimous ratio upstream U alignment.
    d_dN_dS : String/Float
        Non-synonimous/synonimous ratio downstream U alignment.
    changes_dN_dS_u : Int
        Number of CDS mutations upstream U
    changes_dN_dS_d : Int
        Number of CDS mutations downstream U
        
    '''
    # before measuring dN/dS we need to have the cds sequences with the same 
    # size, so we add the gaps.
    subj_cds, query_cds = make_aligned_cds(row['Subj_CDS'], row['Query_CDS'], 
                                           row['Subj_align_prot_seq'], 
                                           row['Q_align_prot_seq'])
    
    u_subj = row['Subj_align_prot_seq'].index('U')
    u_query = row['Q_align_prot_seq'].index('U')
    # calculate the index of the good U, in case there are two or more 
    # aligned U's.
    u_subj, u_query = UGA(row['Subj_align_prot_seq'], 
                          row['Q_align_prot_seq'],
                          u_subj, u_query)
    # 
    changes_dN_dS_u = count_coding_changes(query_cds[:u_query * 3], 
                                           subj_cds[:u_subj * 3])
    changes_dN_dS_d = count_coding_changes(query_cds[(u_query + 1) * 3:], 
                                           subj_cds[(u_subj + 1) * 3:])
    # 
    tupla_dN_dS_u = count_coding_sites(query_cds[:u_query * 3])
    tupla_dN_dS_d = count_coding_sites(query_cds[(u_query + 1) * 3:])
    
    if tupla_dN_dS_u[0] == 0:
        udN = 'NA'
    else:
        udN = changes_dN_dS_u[0]/tupla_dN_dS_u[0]
    if tupla_dN_dS_u[1] == 0:
        udS = 'NA'
    else:
        udS = changes_dN_dS_u[1]/tupla_dN_dS_u[1]
    
    if tupla_dN_dS_d[0] == 0:
        ddN = 'NA'
    else:
        ddN = changes_dN_dS_d[0]/tupla_dN_dS_d[0]
    if tupla_dN_dS_d[1] == 0:
        ddN = 'NA'
    else:
        ddS = changes_dN_dS_d[1]/tupla_dN_dS_d[1]
    
    if udN == 0 or udN == 'NA' or udS == 0 or udS == 'NA':
        u_dN_dS = 'NA'
    else:
        u_dN_dS = round(udN/udS, 4)
    if ddN == 0 or ddN == 'NA' or ddS == 0 or ddS == 'NA':
        d_dN_dS = 'NA'
    else:
        d_dN_dS = round(ddN/ddS, 4)
    
    return u_dN_dS, d_dN_dS, changes_dN_dS_u, changes_dN_dS_d

def alignments_blastp(selenocandidates, db, n_cpu, 
                      blastp_outfile, fasta_query_outfile):
    '''
    This function runs blastp, it is used to get the name of the proteins 
    that have passed all the filters.
    
    Parameters
    ----------
    selenocandidates : Dataframe
        Dataframe with the selenoprotein candidates.
    db : Uniprot database, easyterm object
        opt['uniprot']
    n_cpu : Easyterm object
        opt['c']
    blastp_outfile : String
        Path for the outfile of the blastp.
    fasta_query_outfile : String
        Path for the fasta file.
        
    Raises
    ------
    Exception
        We raise an exception to stop the program in case returncode returns 
        different than zero, indicating that subprocess.run hasn't run 
        successfully.
        We also print the stdout and stderr to know more about the problem.
        
    Returns
    -------
    blastp_outfile : String
        Path for the outfile of the blastp
        
    '''
    
    write(f'Running blastp, results at {blastp_outfile}')
    
    fasta_outfile = ''
    
    with open(fasta_query_outfile, 'w') as fw:
        for i, row in selenocandidates.iterrows():
            q_no_hypens = row['Q_align_prot_seq'].replace('-', '')
            fw.write(f">{row['Q_ID']}\n")
            fw.write(f"{q_no_hypens}\n")
    # command to run tblastx with the query, subject and outfile
    format_6_cmd = (
        'blastp -task blastp-fast -query ' + fasta_query_outfile + 
        ' -db ' + db + ' -out ' + blastp_outfile + 
        ' -num_threads ' + str(n_cpu) + ' -max_hsps ' + str(10) + 
        " -outfmt '6 sacc stitle sstart send sframe qacc " + 
        "qstart qend qframe bitscore evalue sseq qseq'")
    # this function shlex.split, split the string into a shell-like syntax
    format_6_cmd_list = shlex.split(format_6_cmd)
    # subprocess module allows you to spawn new processes, connect to their 
    # input/output/error pies, and obtain their return codes.
    x = subprocess.run(format_6_cmd_list, capture_output=True) 
    
    if x.returncode != 0:
        print(x.stderr, x.stdout)
        raise Exception()
        
    remove(fasta_query_outfile)
    
    # creates a table DataFrame and names the columns
    uniprot_IDs = pd.read_table(blastp_outfile, 
                                names=['Subj_ID', 'Uniprot_Title', 
                                       'Subj_align_s', 'Subj_align_e', 
                                       'Subj_fr', 'Q_ID', 'Q_align_s',
                                       'Q_align_e', 'Q_fr', 'Score', 
                                       'Evalue', 'Subj_align_prot_seq', 
                                       'Q_align_prot_seq'])
    # 
    uniprot_IDs = uniprot_IDs.groupby(
        'Q_ID', as_index = False).agg({'Uniprot_Title': join_titles})
    selenoproteins = selenocandidates.join(
        uniprot_IDs.set_index('Q_ID'), on='Q_ID')
    
    return selenoproteins, uniprot_IDs

def join_titles(uniprot_list):
    '''
    Generates the Uniprot titles.
    
    Parameters
    ----------
    uniprot_list : String
        List of the different titles of the sequences annotated in the uniprot
        database.

    Returns
    -------
    out : String
        All the different titles (max. 10) that the sequence have in the 
        uniprot database.
    
    '''
    
    uniprot_list = uniprot_list[:10]
    already_ind = set()
    out = ''
    
    for title in uniprot_list:
        short_title = ''
        words = title.split()
        for word in words:
            if word.startswith('UniRef'):
                continue
            elif word.startswith('n='):
                break
            else:
                short_title += word + ' '
        if short_title not in already_ind:
            out += short_title.strip() + ', '
            already_ind.add(short_title)
        
    return out[:-2]

def evo_conservation(selenoproteins, dictionary_sel):
    '''
    Third and last filter of the script.
    
    This function is used to measure the evolutionary conservation between 
    both query and subject protein sequences.
    
    Parameters
    ----------
    selenoproteins_candidates : Dataframe
        Dataframe with all the selenoprotein candidates.
    dictionary_sel : List of dictionaries
        List of dictionaries with the values of the BLOSUM Matrix.
        
    Returns
    -------
    outfile : String
        Outfile for the comparison file (blast default format).
    selenoproteins : Dataframe
        Final dataframe with the selenoprotein candidates.
        
    '''
    
    write(f'Evo-conservation of the sequences')
    
    outfile = ''
    
    list_up_dN_dS = list()
    list_down_dN_dS = list()
    
    selenocandidates = pd.DataFrame()
    
    for i, row in selenoproteins.iterrows():
        
        (u_dN_dS_grade, d_dN_dS_grade, 
         u_dN_dS_changes, d_dN_dS_changes) = dN_dS(row)
        up_changes = u_dN_dS_changes[0] + u_dN_dS_changes[1]
        down_changes = d_dN_dS_changes[0] + d_dN_dS_changes[1]
        if up_changes <= 5 or down_changes <= 5:
            continue
        
        list_up_dN_dS.append(u_dN_dS_grade)
        list_down_dN_dS.append(d_dN_dS_grade)
        
        selenocandidates = selenocandidates.append(selenoproteins.iloc[[i]], 
                                                   ignore_index=True)
        comparison_string = ''
        Us_string = ''
        
        for index, x in enumerate(row['Q_align_prot_seq']):
            if x == '-' or row['Subj_align_prot_seq'][index] == '-':
                comparison_string += ' '
                Us_string += ' '
                continue
            
            if x == row['Subj_align_prot_seq'][index]:
                comparison_string += x
            elif dictionary_sel[x][row['Subj_align_prot_seq'][index]] > 0:
                comparison_string += '+'
            else:
                comparison_string += ' '
                
            if x == 'U' or row['Subj_align_prot_seq'][index] == 'U':
                Us_string += '^'
            else:
                Us_string += ' '
        
        n = 60
        comparison_chunks = [comparison_string[y:y+n] 
                             for y in range(0, len(comparison_string), n)]
        query_chunks = [row['Q_align_prot_seq'][q:q+n] 
                        for q in range(0, len(comparison_string), n)]
        subj_chunks = [row['Subj_align_prot_seq'][s:s+n] 
                       for s in range(0, len(comparison_string), n)]
        U_chunks = [Us_string[u:u+n] 
                       for u in range(0, len(comparison_string), n)]
        
        outfile += f'\n'
        outfile += f" ID = {row['Attribute']},  Score = {row['Score']},  Evalue = {row['Evalue']},  Density Score = {row['Density Score']}\n"
        outfile += f" Subj_ID = {row['Chromosome']},  Query_ID = {row['Q_ID']}\n"
        outfile += f" Uniprot_ID = {row['Uniprot_Title']}\n"
        outfile += f" udN/dS = {u_dN_dS_grade}, ddN/dS = {d_dN_dS_grade}\n"
        outfile += f" uNSC + uSC = {up_changes} , dNSC + dSC = {down_changes}\n"
        outfile += f'\n'
        
        gaps_query = 0
        gaps_subj = 0
        for idx, chunk in enumerate(comparison_chunks):
            gaps_query += query_chunks[idx].count('-')
            gaps_subj += subj_chunks[idx].count('-')
            
            if idx == 0:
                if row['Q_Strand'] == '+':
                    outfile += f"Query  {row['Q_align_s'] + n * idx * 3:<5d}  {query_chunks[idx]}  {row['Q_align_s'] + (n * idx + len(query_chunks[idx]) - gaps_query) * 3:<5d}\n"
                else:
                    outfile += f"Query  {row['Q_align_e'] - n * idx * 3:<5d}  {query_chunks[idx]}  {row['Q_align_e'] + (-1 * (n * idx + len(query_chunks[idx])) + gaps_query) * 3:<5d}\n"
                acumulated_gaps_query = gaps_query
            else:
                if row['Q_Strand'] == '+':
                    outfile += f"Query  {row['Q_align_s'] + n * idx * 3 + 1 - acumulated_gaps_query * 3:<5d}  {query_chunks[idx]}  {row['Q_align_s'] + (n * idx + len(query_chunks[idx]) - gaps_query) * 3:<5d}\n"
                else:
                    outfile += f"Query  {row['Q_align_e'] - n * idx * 3 - 1 + acumulated_gaps_query * 3:<5d}  {query_chunks[idx]}  {row['Q_align_e'] + (-1 * (n * idx + len(query_chunks[idx])) + gaps_query) * 3:<5d}\n"
                acumulated_gaps_query = gaps_query
                
            outfile += f'              {chunk}\n'
            
            if idx == 0:
                if row['Strand'] == '+':
                    outfile += f"Sbjct  {row['Start'] + n * idx * 3:<5d}  {subj_chunks[idx]}  {row['Start'] + (n * idx + len(subj_chunks[idx]) - gaps_subj) * 3:<5d}\n"
                else:
                    outfile += f"Sbjct  {row['End'] - n * idx * 3:<5d}  {subj_chunks[idx]}  {row['End'] + (-1 * (n * idx + len(subj_chunks[idx])) + gaps_subj) * 3:<5d}\n"
                acumulated_gaps_subj = gaps_subj
            else:
                if row['Strand'] == '+':
                    outfile += f"Sbjct  {row['Start'] + n * idx * 3 + 1 - acumulated_gaps_subj * 3:<5d}  {subj_chunks[idx]}  {row['Start'] + (n * idx + len(subj_chunks[idx]) - gaps_subj) * 3:<5d}\n"
                else:
                    outfile += f"Sbjct  {row['End'] - n * idx * 3 - 1 + acumulated_gaps_subj * 3:<5d}  {subj_chunks[idx]}  {row['End'] + (-1 * (n * idx + len(subj_chunks[idx])) + gaps_subj) * 3:<5d}\n"
                acumulated_gaps_subj = gaps_subj
                
            outfile += f'              {U_chunks[idx]}\n'
        
    selenocandidates['dN/dS upstream'] = list_up_dN_dS
    selenocandidates['dN/dS downstream'] = list_down_dN_dS
    
    return outfile, selenocandidates

def outputs(candidates, cds_q, path_cds_q, cds_t, path_cds_t, pep_q, 
            path_pep_q, pep_t, path_pep_t, gff_q, gff_t, path_query_gff, 
            path_subj_gff, dotplot, path_dotplot):
    '''
    Function with the several optional script outfiles dependent on boolean 
    values.
    
    Parameters
    ----------
    candidates : Dataframe
        Dataframe with all the selenoprotein candidates.
    cds_q : Boolean, easyterm object
        opt['cds_q']
    path_cds_q : String
        Path to save the cds sequence from the query if opt['cds_q'] == True.
    cds_t : Boolean, easyterm object
        opt['cds_t']
    path_cds_t : String
        Path to save the cds sequence from the target if opt['cds_t'] == True.
    pep_q : Boolean, easyterm object
        opt['pep_q']
    path_pep_q : String
        Path to save the protein sequence from the query if opt['pep_q'] == True.
    pep_t : Boolean, easyterm object
        opt['pep_t']
    path_pep_t : String
        Path to save the protein sequence from the target if opt['pep_t'] == True.
    gff_q : Boolean, easyterm object
        opt['dff_q']
    gff_t : Boolean, easyterm object
        opt['dff_t']
    path_query_gff : String
        Path to save the gff file from the query if opt['dff_q'] == True.
    path_subj_gff : String
        Path to save the gff file from the target if opt['dff_t'] == True.
    dotplot : Boolean, easyterm object
        opt['dotplot']
    path_dotplot : String
        Path to save the dotplot if opt['dotplot'] == True.
        
    Returns
    -------
    None
    
    '''
    write(f'Producing output')
    
    output_cds_q = ''
    output_cds_t = ''
    output_pep_q = ''
    output_pep_t = ''
    
    for i, row in candidates.iterrows():
        subj_cds, query_cds = make_aligned_cds(row['Subj_CDS'], 
                                               row['Query_CDS'], 
                                               row['Subj_align_prot_seq'], 
                                               row['Q_align_prot_seq'])
        if cds_q == True:
            output_cds_q += f"<{row['Q_ID']}\n"
            output_cds_q += f"{query_cds}\n"
        if cds_t == True:
            output_cds_t += f"<{row['Chromosome']}\n"
            output_cds_t += f"{subj_cds}\n"
        if pep_q == True:
            output_pep_q += f"<{row['Q_ID']}\n"
            output_pep_q += f"{row['Q_align_prot_seq']}\n"
        if pep_t == True:
            output_pep_t += f"<{row['Chromosome']}\n"
            output_pep_t += f"{row['Subj_align_prot_seq']}\n"
    
    if cds_q == True:
        with open(path_cds_q, 'w') as fw:
            fw.write(output_cds_q)
    if cds_t == True:
        with open(path_cds_t, 'w') as fw:
            fw.write(output_cds_t)
    if pep_q == True:
        with open(path_pep_q, 'w') as fw:
            fw.write(output_pep_q)
    if pep_t == True:
        with open(path_pep_t, 'w') as fw:
            fw.write(output_pep_t)
    
    if dotplot == True:
        dot_plot(candidates, path_dotplot)
    
    write(f'Dividing columns into query and subject dataframes')
    
    candidates.rename(columns = {'Attribute':'Gene_ID'}, inplace=True)
    candidates['Feature'] = 'CDS'
    candidates['Source'] = 'Twinstop'
    
    if gff_q == True:
        query_df = candidates.copy()
        query_df.drop(['Chromosome', 'Start', 'End', 'Strand'], 
                      axis=1, inplace=True)
        query_df = query_df.rename(
            columns = {'Q_ID':'Chromosome', 'Q_align_s':'Start',
                       'Q_align_e':'End', 'Q_Strand':'Strand'})
        query_df_reduced = query_df.loc[
            :, ['Chromosome', 'Source', 'Feature', 'Start', 'End', 
                'Strand', 'Score', 'Gene_ID', 'Uniprot_Title']]
            
        py_query = pr.PyRanges(query_df_reduced)
        py_query.to_gff3(path = path_query_gff)
        
    if gff_t == True:
        candidates_reduced = candidates.loc[
            :, ['Chromosome', 'Source', 'Feature', 'Start', 'End', 
                'Strand', 'Score', 'Gene_ID', 'Uniprot_Title']]
            
        py_subj = pr.PyRanges(candidates_reduced)
        py_subj.to_gff3(path = path_subj_gff)
        
    write(f'Done')

def main():
    
    tracemalloc.start()
    now = datetime.now()
    current_time = now.strftime('%H:%M:%S')
    
    __version__ = '0.0.1'
    
    help_msg = """This program allows us to make inputs to the pipeline
    directly from the terminal.
    
    ### Input/Output:
    -q : query file
    -t : db file
    -o : output path
    -c : number of CPUs used to run BLAST
    -b : control the run of default tBLASTx
    -f : force the run of tBLASTx in format 6
    -cds_q : control the creation of the nucleotidic sequences of the 
             candidates in the query.
    -cds_t : control the creation of the nucleotidic sequences of the 
                candidates in the subject.
    -dotplot : control the creation of one or more dotplots with the 
               candidates sequences.
    -pep_q : control the creation of the proteic sequences of the 
             candidates in the query.
    -pep_t : control the creation of the nucleotidic sequences of the 
                candidates in the subject.
    -dff_q : control the creation of a .dff file with the annotations of the 
             candidates in the query.
    -dff_t : control the creation of a .dff file with the annotations of the 
                candidates in the subject.
    -uniprot : uniprot blast database
    -cons_up : conservation minimum value upstream
    -cons_down : conservation minimum value downstream
    -n_section : control the rerun of the sections of the scripts
    
    ### Options:
    -print_opt: print currently active options
    -h | --help: print this help and exit"""
    
    def_opt= {
        'q': 'query',
        't': 'subject',
        'o': '/home/ssanchez/seleno_prediction/outputs/',
        'c': 4,
        'd': False,
        'f': False,
        'cds_q': False,
        'cds_t': False,
        'dotplot': False,
        'pep_q': False,
        'pep_t': False,
        'gff_q': False,
        'gff_t': False,
        'uniprot': 
            '/home/ssanchez/seleno_prediction/enlaces_fungi/uniref50.fasta',
        'cons_up': 50,
        'cons_down': 150,
        'n_section': 1000}
        
    opt = command_line_options(def_opt, help_msg)
    
    check_file_presence(opt['q'], 'query')
    # return the base name of pathname path
    query = path.basename(opt['q']).split('.')[0]
    subject = path.basename(opt['t']).split('.')[0]
    
    if path.exists(opt['o']) == False:
        makedirs(opt['o'])
    # Main paths
    resultados = opt['o'] + 'tblastx.tsv'
    path_fragments = opt['o'] + 'all_orfs.tsv'
    selected_IDs = opt['o'] + 'nov_orfs.tsv'
    path_table_df = opt['o'] + 'ext_orfs.txt'
    path_pairwise = opt['o'] + 'aln_orfs.tsv'
    # Temporal paths
    query_gff = opt['o'] + 'table_query.gff'
    subject_gff = opt['o'] + 'table_subj.gff'
    out_subj_gff = opt['o'] + 'out_subj.gff'
    out_query_gff = opt['o'] + 'out_query.gff'
    path_fasta_outfile = opt['o'] + 'fasta_seq.txt'
    # Output paths
    path_selenocandidates = opt['o'] + 'candidates.tsv'
    out_blastp = opt['o'] + 'candidates_blastp.tsv'
    comparison_outfile = opt['o'] + 'candidates_pretty.txt'
    path_cds_q = opt['o'] + 'candidates_query.cds.fa'
    path_cds_t = opt['o'] + 'candidates_target.cds.fa'
    path_pep_q = opt['o'] + 'candidates_query.pep.fa'
    path_pep_t = opt['o'] + 'candidates_target.pep.fa'
    path_gff_q = opt['o'] + 'candidates_query.gff'
    path_gff_t = opt['o'] + 'candidates_target.gff'
    path_dot_plot = opt['o'] + 'candidates_dotplot.png'
    fragments_dot_plot = opt['o'] + 'fragments_dotplot.png'
    overlapping_dot_plot = opt['o'] + 'overlapping_dotplot.png'
    
    write(f'TwinStop {__version__}')
    write(f'{current_time}')
    write(f'\n')
    write(f'### PHASE 1: TBLASTX')
    
    tblastx_hits = run_tblastx(opt['q'], opt['t'], opt['o'], 
                               opt['c'], opt['d'], opt['f'], 
                               query, subject, resultados)
    
    if path.exists(resultados) == False or opt['n_section'] < 2:
        tblastx = conv_pyran(tblastx_hits, opt['t'], opt['q'])
        # when Dataframe, 'path_or_buf' is the argument
        # if PyRanges, it is 'path'.
        tblastx.to_csv(sep='\t', path_or_buf = resultados)
    else:
        write(f'Reading {resultados}')
        tblastx = pd.read_csv(resultados, sep='\t', header=0, 
                              index_col=0)
        if len(tblastx) == 0:
            write(f'Empty file {final_table}')
            
    dictionary_matrix = dictionary_seleno()
    
    write(f'\n')
    write(f'### PHASE 2: FRAGMENTATION')
    
    if path.exists(path_fragments) == False or opt['n_section'] < 3:
        fragments = max_score_block(tblastx, dictionary_matrix)
        del(tblastx)
        fragments.to_csv(sep='\t', path_or_buf = path_fragments)
    else:
        write(f'Reading {path_fragments}')
        fragments = pd.read_csv(path_fragments, sep='\t', header=0, 
                                index_col=0)
        if len(fragments) == 0:
            write(f'Empty file {fragments}')
            
    if opt['dotplot'] == True:
        dot_plot(fragments, fragments_dot_plot)
    
    write(f'\n')
    write(f'### PHASE 3: OVERLAP FILTER')
    
    if path.exists(selected_IDs) == False or opt['n_section'] < 4:
        clusters = overlapping(fragments)
        del(fragments)
        clusters.to_csv(sep='\t', path_or_buf = selected_IDs)
    else:
        write(f'Reading {selected_IDs}')
        clusters = pd.read_csv(selected_IDs, sep='\t', header=0, 
                               index_col=0)
        if len(clusters) == 0:
            write(f'Empty file {clusters}')
            
    if opt['dotplot'] == True:
        dot_plot(clusters, overlapping_dot_plot)
    
    write(f'\n')
    write(f'### PHASE 4: EXTEND ORFS')
    
    if path.exists(path_table_df) == False or opt['n_section'] < 5:
        table_df = extend_orfs(clusters, subject_gff, query_gff, 
                               opt['t'], opt['q'], out_subj_gff, 
                               out_query_gff, opt['o'])
        table_df.to_csv(sep='\t', path_or_buf = path_table_df)
    else:
        write(f'Reading {path_table_df}')
        table_df = pd.read_csv(path_table_df, sep='\t', header=0, 
                               index_col=0)
        if len(table_df) == 0:
            write(f'Empty file {table_df}')
            
    write(f'\n')
    write(f'### PHASE 5: ALIGNMENTS')
    
    if path.exists(path_pairwise) == False or opt['n_section'] < 6:
        joined_df = pairwise(table_df, dictionary_matrix)
        joined_df.to_csv(sep='\t', 
                         path_or_buf = path_pairwise)
    else:
        write(f'Reading {path_pairwise}')
        joined_df = pd.read_csv(path_pairwise, sep='\t', header=0, 
                                index_col=0)
        if len(joined_df) == 0:
            write(f'Empty file {joined_df}')
    
    write(f'\n')
    write(f'### PHASE 6: FILTER')
    
    if path.exists(path_selenocandidates) == False or opt['n_section'] < 7:
        candidates = UGA_alignments(joined_df, dictionary_matrix, 
                                    opt['cons_up'], opt['cons_down'])
        candidates.to_csv(sep='\t', 
                          path_or_buf = path_selenocandidates)
    else:
        write(f'Reading {path_selenocandidates}')
        candidates = pd.read_csv(path_selenocandidates, 
                                 sep='\t', header=0, index_col=0)
        if len(candidates) == 0:
            write(f'Empty file {selenocandidates}')
    
    write(f'\n')
    write(f'### PHASE 7: BLASTP FOR TITLES')
    
    if not 'Uniprot_Title' in candidates.columns:
        candidates, candidates_blastp = alignments_blastp(
            candidates, opt['uniprot'], opt['c'], 
            out_blastp, path_fasta_outfile)
    
    write(f'\n')
    write(f'### PHASE 8: OUTPUTS')
    
    candidates_pretty, candidates = evo_conservation(candidates, 
                                                     dictionary_matrix)
    candidates.to_csv(sep='\t', path_or_buf = path_selenocandidates)
    
    with open(comparison_outfile, 'w') as fw:
        fw.write(candidates_pretty)
    
    outputs(
        candidates, opt['cds_q'], path_cds_q, opt['cds_t'], path_cds_t,
        opt['pep_q'], path_pep_q, opt['pep_t'], path_pep_t, opt['gff_q'], 
        opt['gff_t'], path_gff_q, path_gff_t, opt['dotplot'], path_dot_plot)
    
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()
    
if __name__ == '__main__':
    main()
   