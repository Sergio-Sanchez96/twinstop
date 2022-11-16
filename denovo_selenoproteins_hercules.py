#!/users-d3/EGB-invitado4/miniconda3/envs/twinstop/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 11:55:40 2022

@author: Sergio Sánchez Moragues
"""

import sys
sys.path.insert(0, '/users-d3/EGB-invitado4/software/pyranges/pyranges/pyranges.py')
import os
import numpy as np
import pyranges as pr
import subprocess
import shlex
import pandas as pd
import tracemalloc
import gc #???
from datetime import datetime
from easyterm import command_line_options, check_file_presence, write
from easybioinfo import count_coding_changes, count_coding_sites, translate
from Bio import pairwise2
from block_selection import block_dict, score
from dotplot_twinstop import dot_plot
from chunk_files import iterate_by_chunks

def run_tblastx(q_file, db_file, outfolder, n_cpu, IsDefaultFileCreated,
                IsFormat6FileForced, format_6_outfile):
    '''
    Runs tblastx in default (optionally) and tabular 6 format.
    Columns selected: subject accession(sacc), subject alignment start (sstart), 
    subject alignment end (send), subject frame (sframe), query accession (qacc), 
    query alignment start (qstart), query alignment end (qend), query frame (qframe)
    bitscore (bitscore), evalue (evalue), query aligned protein sequence (qseq), 
    subject aligned protein sequence (sseq).

    Parameters
    ----------
    q_file : String
        Path where the transcriptome from the query is located
    db_file : String
        Path where the transcriptome from the subject is located. It must be 
        previously recognized by Blast with the command 'makeblastdb'.
    outfolder : String, easyterm object
        Path to the output folder (opt['o']).
    n_cpu : Int, easyterm object
        opt['c']
    IsDefaultFileCreated : Boolean value, easyterm object
        opt['b']
    IsFormat6FileForced : Boolean value, easyterm object
        opt['f']
    tblastxF6_output : String
        Path where the tblastx tabular format 6 output will be located

    Raises
    ------
    Exception
        We raise an exception to stop the program in case returncode is
        different that zero, indicating that subprocess.run hasn't run
        successfully.
        We also print the stdout and stderr to know more about the problem

    Returns
    -------
    tblastxF6_output : Tblastx outfile
        With the results of the tblastx hits in the tabular format 6
    '''

    # if IsDefaultFileCreated: # decide about the creation of the default file
        # write(f'Running default tblastx')
        # temp_default_outfile = outfolder + 'temp_default_tblastx_hits.txt'
        # name and location of the default file
        # default_outfile = outfolder + 'default_tblastx_hits.txt'
        # string to run tblastx with the query, subject, outfile and number of
        # cpu used.
        # default_cmd = ('tblastx -q_file ' + q_file + ' -db ' + db_file + ' -out ' +
        #                temp_default_outfile + ' -num_threads ' + str(n_cpu))
        # shlex.split, splits the string into a shell-like syntax
        # default_cmd_list = shlex.split(default_cmd)
        # subprocess module allows you to spawn new processes, connect to 
        # their input/output/error pies, and obtain their return codes.
        # x = subprocess.run(default_cmd_list, capture_output=True)
        # if x.returncode != 0:
        #     print(x.stderr, x.stdout)
        #     raise Exception()
        # exit status of the process. Typically, an exit status of 0 indicates
        # that it ran successfully.
        # in case the program fails, we print the standard output (stdout) and
        # the standard error (stderr) and an exception is raised to stop the program.

        # replace the temp file to be sure than the tblastx is completed
        # os.replace(temp_default_outfile, default_outfile)
        # write(f'Default tblastx ran successfully, results in {default_outfile}')

    # we run the tblastx format 6 table only if it does not exist or if we 
    # force it.
    if not os.path.exists(format_6_outfile) or IsFormat6FileForced:
        write(f'Running format 6 tblastx')
        temp_format_6_outfile = outfolder + 'temp_tblastx.tsv'
        # command to run tblastx with the specific columns we want
        format_6_cmd = (
            'tblastx -query ' + q_file + ' -db ' + db_file +
            " -outfmt '6 sacc sstart send sframe qacc qstart" + 
            " qend qframe bitscore evalue sseq qseq' -out " + 
            temp_format_6_outfile + ' -num_threads ' + str(n_cpu))
        # print(sys.getsizeof(format_6_cmd))
        format_6_cmd_list = shlex.split(format_6_cmd)
        # print(sys.getsizeof(format_6_cmd_list))
        y = subprocess.run(format_6_cmd_list, capture_output=True)
        # print(sys.getsizeof(y))
        if y.returncode != 0:
            print(y.stderr, y.stdout)
            raise Exception()

        os.replace(temp_format_6_outfile, format_6_outfile)
        write(f'Format 6 tblastx ran successfully, results in {format_6_outfile}')

def pandas1(tblastx_outfile, tblastx_df_file):
    '''
    Converts the tblastx output in a pandas dataframe and divides it into
    chunks to reduce memory usage, then they are joined back.

    Parameters
    ----------
    tblastx_outfile : String
        Format 6 output file with the resulting tblastx hits
    n_chunks : Int, easyterm object
        opt['n_chunks']
    db_file : String
        Path where the transcriptome from the subject is located.
    q_file : String
        Path where the transcriptome from the query is located

    Returns
    -------
    None
    '''

    write(f'Converting {tblastx_outfile} into a pandas DataFrame')
    # creates a DataFrame and names the columns
    colnames = ['Chromosome', 'Start', 'End', 'Subj_fr', 'Q_ID', 'Q_align_s', 'Q_align_e',
               'Q_fr', 'Score', 'Evalue', 'Subj_align_prot_seq', 'Q_align_prot_seq']
    df = pd.read_table(tblastx_outfile, names=colnames)
    # range() does not count the last nº, that is why we put a +1
    df['ID'] = range(
        1, len(df) + 1) # Identification column
    write(f'ID column created')
    colnames.append('ID')
    df.to_csv(sep='\t', path_or_buf=tblastx_df_file)
    # df.to_csv(sep='\t', path_or_buf=tblastx_df_file, index=False)
    # Default: header=True, index=True
    # table_subj.to_csv(sep='\t', path_or_buf = tblastx_df, index = False)
    # if opt['chunksize'] == 0:
    #     chunksize = table_subj.shape[0]//opt['n_chunks']
    # else:
    #     chunksize = opt['chunksize']
    return colnames

def chunking(fname, n_chunks, func, colnames, chunks_file):
    '''
    '''
    # np.array_split() divides the df in n number of chunks (returns a list)
    # for i, chunk in enumerate(np.array_split(df, n_chunks)):
    iterator = iterate_by_chunks(fname, nchunks=n_chunks)
    # iterator = pd.read_csv(tblastx_df, sep='\t',chunksize=chunksize, header=0)
    # for i, chunk in enumerate(iterator):
    for chunkindex in range(n_chunks):
        if chunkindex == 0:
            header = 0
        else:
            header = None
        chunkdf = pd.read_csv(iterator, engine='python', header=header, names=colnames, sep='\t')
        df = func(chunkdf)
        # if chunkindex == 0:
        #     mode = 'w'
        #     header = chunkindex
        # else:
        #     mode = 'a'
        mode = 'w' if chunkindex == 0 else 'a'
        header = chunkindex == 0
        # when Dataframe, 'path_or_buf' is the argument | if PyRanges, it is 'path'
        df.to_csv(sep='\t', path_or_buf=chunks_file, header=header, mode=mode, index=False)
        del df # important to delete df variable before starting next loop

def conv_pandas_df(table_subj):
    '''
    Divides the tblastx dataframe in two pandas DataFrames, one for subject
    the other for query.

    Parameters
    ----------
    table_subj : Dataframe
        Tblastx tabular format 6 dataframe

    Returns
    -------
    table_query : Dataframe
        Dataframe only with the query-related columns.
    table_subj : Dataframe
        Dataframe only with the subject-related columns.
    '''

    # where 'Subj_fr' is greater than 0 it will be put a '+' in the 'Strand' column,
    # if not, a '-' (Pyranges format).
    table_subj['Strand'] = (
        (table_subj['Subj_fr'] > 0).replace(
            {True: '+', False: '-'}))
    write(f'Strand column created')
    # creates a boolean Series
    indexer = table_subj['Start'] > table_subj['End']
    indexer2 = table_subj['Q_align_s'] > table_subj['Q_align_e']
    # switches the columns of 'Start' and 'End' of both the query and the 
    # subject if they are in a negative frame.
    table_subj.loc[
        indexer, ['Start', 'End']] = table_subj.loc[
            indexer, ['End', 'Start']].values
    write(f'Indexers done')
    # substitutes the value of the column only where indexer is True.
    table_subj.loc[
        indexer2, ['Q_align_s', 'Q_align_e']] = table_subj.loc[
            indexer2, ['Q_align_e', 'Q_align_s']].values
    # BLAST is a 1-base program (does not take the last value), #???
    # while Pyranges is 0-base (takes the first and last value).
    table_subj['Start'] = table_subj['Start'] - 1
    table_subj['Q_align_s'] = table_subj['Q_align_s'] - 1
    # print(sys.getsizeof(table_subj))
    # print(table_subj.memory_usage(deep=True))
    write(f'Start values adapted to 0-base program')
    # divide the dataframe into two dataframes (one for the subject and the 
    # other for the query), ready to be transformed into Pyranges.

    return query_subject_dfs(table_subj) # we return both dataframes

def query_subject_dfs(table_subj):
    '''
    Divides the whole dataframe into two, one for the query and the 
    other for the subject, with the format needed to convert them into 
    Pyranges objects.

    Parameters
    ----------
    table_subj : Dataframe
        Dataframe with all the columns about the blasthits

    Returns
    -------
    table_query : Dataframe
        Dataframe only with the query-related columns
    table_subj : Dataframe
        Dataframe only with the subject-related columns
    '''

    write(f'Dividing columns into query and subject dataframes')
    table_query = table_subj.copy() # creates table_query from the table_subj
    # drops the query-related columns
    table_subj.drop(['Q_ID', 'Q_align_s', 'Q_align_e', 'Q_fr',
                     'Q_align_prot_seq'], axis=1, inplace=True)
    write(f'Dropping Query columns in subject df')
    # renames the columns to fit into the Pyranges format (Chromosome, Start, End, Strand (+/-))
    table_query = table_query.rename(
        columns={'Chromosome': 'Subj_ID', 'Start': 'Subj_align_s', 
                 'End': 'Subj_align_e', 'Q_ID': 'Chromosome', 
                 'Q_align_s': 'Start', 'Q_align_e': 'End'})
    write(f'Renaming Subject columns in query df')
    # drops the subject-related columns ('Score' and 'Evalue' are left in 
    # the subj_table).
    table_query.drop(['Strand', 'Subj_ID', 'Subj_align_s', 'Subj_align_e', 'Subj_fr',
                      'Subj_align_prot_seq', 'Score', 'Evalue'], axis=1, inplace=True)
    write(f'Dropping Subject columns in query df')
    table_query['Strand'] = table_query['Q_fr'].copy()
    # Strand column needs to have '+' or '-' only
    table_query['Strand'] = (
        table_query['Strand'] > 0).replace({True:'+', False:'-'})
    write(f'Replacing Strand column by +/- in query df')
    # print(sys.getsizeof(table_subj))
    # print(sys.getsizeof(table_query))
    # we return both the query and subject dataframes.
    write(f'Returning query and subject dataframes')

    return table_query, table_subj

def conv_pyran(tblastx_df, db_file, q_file):
    '''
    Converts dataframes into pyranges, extract the CDS and protein sequences 
    from both query and subject, and then rejoins them into one (after 
    transforming into dataframes again).

    Parameters
    ----------
    tblastx_df : Dataframe
        Dataframe with the tblastx hits.
    db_file : String
        Path where the transcriptome from the subject is located.
    q_file : String
        Path where the transcriptome from the query is located.

    Returns
    -------
    final_table_df : Dataframe
        Dataframe table with the tblastx columns plus the results of the 
        pyranges CDS sequences.
    '''

    query_table, subject_table = conv_pandas_df(tblastx_df)
    del tblastx_df
    write(f'Converting subject and query dataframes into pyranges')
    query_table, subject_table = translate_prot(subject_table, query_table, 
                                                db_file, q_file)
    # print(sys.getsizeof(subject_table))
    # print(sys.getsizeof(query_table))
    query_table.drop(['Strand'], axis=1, inplace=True)
    # we need to rename query's columns before joining back the two dataframes
    query_table = query_table.rename(columns={'Chromosome': 'Q_ID', 
                                              'Start': 'Q_align_s', 
                                              'End': 'Q_align_e'})
    # joins subject and query DataFrames according to ID column
    # set_index() drop=True by default
    final_table_df = subject_table.join(query_table.set_index('ID'), on='ID')
    del subject_table
    del query_table
    # final_table_df = pd.merge(subject_table, query_table, left_index=True, 
    #                           right_index=True)
    final_table_df = final_table_df.reindex(
        columns=['ID', 'Chromosome', 'Start', 'End', 'Strand', 'Subj_fr', 
                 'Q_ID','Q_align_s', 'Q_align_e', 'Q_fr', 'Score', 
                 'Evalue', 'Subj_align_prot_seq', 'Q_align_prot_seq'])
    final_table_df.sort_values(by='ID', inplace=True, ignore_index=True)
    # print(sys.getsizeof(final_table_df))
    # print(final_table_df.memory_usage(deep=True))
    return final_table_df # returns the joined dataframe

def max_score_block(selected_IDs, dictionary_matrix):
    '''
    Selection of the fragment (between two stops) with the best score of each
    alignment.

    Parameters
    ----------
    selected_IDs : Dataframe
        Blast hits dataframe
    dictionary_matrix : Dictionary of tuples
        Dictionary with the Matrix BLOSUM62 values

    Returns
    -------
    selected_IDs : Dataframe
        Dataframe updated with the selected fragment from each tblastx-hit
    '''

    # lists to update the general dataframe
    list_subj_start = list()
    list_subj_end = list()
    list_query_start = list()
    list_query_end = list()
    list_query_prot = list()
    list_subj_prot = list()
    list_score = list()

    write(f'Taking the fragments with the best score')

    # iters the dataframe by rows
    for i, row in selected_IDs.iterrows():
        # selection of the fragment with the highest score
        max_score_frag = block_dict(row['Q_align_prot_seq'], 
                                    row['Subj_align_prot_seq'], 
                                    dictionary_matrix)

        list_query_prot.append(row['Q_align_prot_seq'][max_score_frag['Align_Start']: max_score_frag['Align_End']])
        list_subj_prot.append(row['Subj_align_prot_seq'][max_score_frag['Align_Start']: max_score_frag['Align_End']])
        list_score.append(max_score_frag['Score'])
        # for the positive frames
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
        # gc.collect()

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
    Creates a dictionary of tuples with the values of the BLOSUM62 Matrix

    Returns
    -------
    dictionary_sel : Dictionary
        Dictionary of tuples with the values of the BLOSUM62 Matrix
    '''

    write(f'Writing the dictionary BLOSUM62')
    dictionary_sel = dict()
    with open('/users-d3/EGB-invitado4/seleno_prediction/data/Matrix_BLOSUM62sel.txt', 'r') as fr:
        for index, row in enumerate(fr):
            # creates a list, using ' ' as sep
            spt = row.split(' ')
            # deletes blank spaces
            spt = list(filter(None, spt))
            if index == 0:
                # delete empthy spaces
                header = [x.strip() for x in spt]
                continue
            # deletes '\n' characters
            spt = [x.strip() for x in spt]
                # converts the values (string) into integers
            ints = [int(x) for x in spt[1:]]
            keys = [(spt[0], aa) for aa in header]

            for ik, k in enumerate(keys):
                dictionary_sel[k] = ints[ik]
            
            # gc.collect()
    return dictionary_sel

def overlapping(real_final_table, selected_IDs):
    '''
    First filter of the script. Based on selecting only the best score among 
    the overlapping hits.

    Parameters
    ----------
    real_final_table : Dataframe
        With all the columns
    selected_IDs : String
        Path where we dataframe will be located after the overlapping filter

    Returns
    -------
    None
    '''

    write(f'Overlapping')
    # converts into PyRanges
    final_table_pr = pr.PyRanges(real_final_table)
    del real_final_table
    # creates a 'Cluster' column identifying the overlapping sequences
    final_table_pr = final_table_pr.cluster(strand=False, slack=0)
    # converts into Dataframe
    final_table_df = final_table_pr.as_df()
    del final_table_pr
    # discards the alignments with evalues greater than 0.05
    final_table_df = final_table_df[final_table_df['Evalue'] < 0.05]
    # creates a series with as many rows as clusters in the dataframe
    cluster_ID = final_table_df['Cluster'].unique()

    for i in cluster_ID:
        # df with only those rows where 'Cluster' column is equal to 'i'
        cluster = final_table_df[final_table_df['Cluster'] == i]
        # keeps the row with the best score
        cluster = cluster[cluster['Score'] == cluster['Score'].max()]
        # substitutes all the cluster by the selected row
        final_table_df[final_table_df['Cluster'] == i] = cluster
        del cluster
        # selected_IDs = selected_IDs.append(cluster, ignore_index=True)
        # gc.collect()
    # replaces the white spaces by NaN
    final_table_df.replace('', np.nan, inplace=True)
    # drops the NaN rows 
    final_table_df.dropna(inplace=True)
    # changes the types of the columns, NaN is only available in float category
    final_table_df = final_table_df.astype({'ID': np.int32, 'Start': np.int32, 'End': np.int32, 'Subj_fr': np.int32,
                                            'Q_align_s': np.int32, 'Q_align_e': np.int32, 'Q_fr': np.int32})
    final_table_df.drop('Cluster', axis=1, inplace=True)
    # returns the Dataframe with only the rows with the best scores 
    # among the overlapping hits.
    final_table_df.to_csv(sep='\t', path_or_buf=selected_IDs, index=False)

def conv_gff(clusters, query_gff, subject_gff):
    '''
    Converts dataframes into gffs

    Parameters
    ----------
    clusters : Dataframe
        General Dataframe with all the columns
    query_gff : String
        Path where the gff of the query will be saved
    subject_gff : String
        Path where the gff of the subject will be saved

    Returns
    -------
    None
    '''

    write(f'GFF format')
    table_query, table_subj = query_subject_dfs(clusters)
    # we have to change 'ID' column's name to fit with the gff format
    table_query = table_query.rename(columns={'ID': 'Attribute'})
    # we also change 'Score' column's name to include it in the 'Attribute' 
    # column.
    # 'extend_orfs.py' does not take the score values for the 'Score' column
    table_subj = table_subj.rename(columns={'ID': 'Attribute',
                                            'Score': 'Value'})
    # creates a new column called 'Feature' with the same value ('CDS') 
    # in all rows.
    table_query['Feature'] = 'CDS'
    table_subj['Feature'] = 'CDS'
    # we only take the columns we need to fulfill the gff format
    # 'Attribute' and 'Q_fr' will be put together in the 'Attribute' column
    table_query_reduced = table_query.loc[:, ['Chromosome', 'Start', 'End', 
                                              'Strand', 'Attribute', 
                                              'Feature', 'Q_fr']]
    # 'Attribute', 'Value' and 'Evalue' will be put together in the 
    # 'Attribute' column.
    table_subj_reduced = table_subj.loc[:, ['Chromosome', 'Start', 'End', 
                                            'Strand', 'Attribute', 'Feature', 
                                            'Value', 'Evalue', 'Subj_fr']]

    py_query = pr.PyRanges(table_query_reduced) # converts into PyRanges
    py_subj = pr.PyRanges(table_subj_reduced)

    py_query.to_gff3(path=query_gff) # converts pyranges object into gff
    py_subj.to_gff3(path=subject_gff)

def extend_orfs(clusters, subject_gff, query_gff, subject, query, out_subj_gff, out_query_gff):
    '''
    Runs Marco's script 'extend_orfs.py' which extends the CDS sequences from 
    the query and the subject both downstream and upstream, until a stop is 
    found or the transcript finishes (downstream) or until an initial codon is 
    found (Methionine) or the transcript finishes (upstream).

    Parameters
    ----------
    clusters : Dataframe
        Dataframe with the best scored tblastx hits among the overlapping ones
    subject_gff : String
        Path where the gff file related to the subject is saved
    query_gff : String
        Path where the gff file related to the query is saved
    subject : String
        Path where the transcriptome of the subject is located
    query : String
        Path where the transcriptome of the query is located
    out_subj_gff : String
        Path where the gff file related to the subject will be saved
    out_query_gff : String
        Path where the gff file related to the query will be saved

    Raises
    ------
    Exception
        We raise an exception to stop the program in case returncode returns 
        different from zero, indicating that subprocess.run hasn't run
        successfully.
        We also print the stdout and stderr to know more about the problem.

    Returns
    -------
    joined_df : Dataframe
        Dataframe with the alignments after extending the orfs.
    '''

    write(f'Extending orfs')
    conv_gff(clusters, query_gff, subject_gff)
    # we run 'extend_orfs' per subject and per query
    # the extension stops only when it encounters an 'TAG' or 'TAA' stop
    subj = ('/users-d3/EGB-invitado4/scripts/extend_orfs.py -i ' + subject_gff + 
            ' -t ' + subject + ' -o ' + out_subj_gff + " -stops 'TAG,TAA'")
    subj_list = shlex.split(subj)
    x = subprocess.run(subj_list, capture_output=True)
    if x.returncode != 0:
        print(x.stderr, x.stdout)
        raise Exception()

    q = ('/users-d3/EGB-invitado4/scripts/extend_orfs.py -i ' + query_gff + ' -t ' + 
         query + ' -o ' + out_query_gff + " -stops 'TAG,TAA'")
    q_list = shlex.split(q)
    y = subprocess.run(q_list, capture_output=True)
    if y.returncode != 0:
        print(y.stderr, y.stdout)
        raise Exception()

    os.remove(subject_gff) # removes the temporal files
    os.remove(query_gff)
    # name the columns of the gff resulting file converting it into a dataframe
    subj_df_gff = pd.read_csv(out_subj_gff, sep='\t', names=['Chromosome', 'Source', 'Feature', 'Start', 'End',
                                                             'Score', 'Strand', 'Frame', 'Attribute'])
    # all the columns that do not have one of this names will be put together in the 'Attribute' column
    query_df_gff = pd.read_csv(out_query_gff, sep='\t', names=['Chromosome', 'Source', 'Feature', 'Start', 'End',
                                                               'Score', 'Strand', 'Frame', 'Attribute'])
    os.remove(out_subj_gff)
    os.remove(out_query_gff)
    # we need to subtract 1 again to the 'start' because 'extend_orfs' has added
    # up 1 directly in each row.
    query_df_gff['Start'] = query_df_gff['Start'] - 1
    subj_df_gff['Start'] = subj_df_gff['Start'] - 1

    return pandas2(subj_df_gff, query_df_gff, subject, query)

def translate_prot(subj_df, query_df, subject, query, CDS_sequences=False):
    '''
    Function to translate the nucleotide sequences into protein using 
    translate(), from Marco's easyterm module.

    Parameters
    ----------
    subj_df : Dataframe
        Dataframe with the subject-related columns
    query_df : Dataframe
        Dataframe with the query-related columns
    subject : String
        Path to the file of the subject
    query : String
        Path to the file of the query

    Returns
    -------
    query_df : Dataframe
        Dataframe  with the query-related columns after translation
    subj_df : Dataframe
        Dataframe with the subject-related columns after translation
    '''

    write(f'Translating into protein')
    query_pr = pr.PyRanges(query_df) # converts into PyRanges
    subj_pr = pr.PyRanges(subj_df)
    del query_df
    del subj_df
    # gc.collect() # deletes del() objects
    if CDS_sequences:
        # gets the CDS sequences
        query_pr.Query_CDS = pr.get_fasta(query_pr, query)
        query_pr.Q_align_prot_seq = [translate(s, genetic_code='1+U') for s in query_pr.Query_CDS]
        subj_pr.Subj_CDS = pr.get_fasta(subj_pr, subject)
        subj_pr.Subj_align_prot_seq = [translate(s, genetic_code='1+U') for s in subj_pr.Subj_CDS]
        write(f'CDS sequences saved')
    else:
        # translates the CDS sequences into protein (conserves the 'U's)
        query_pr.Q_align_prot_seq = (
            [translate(s, genetic_code='1+U') for s in pr.get_fasta(query_pr, query)])
        subj_pr.Subj_align_prot_seq = (
            [translate(s, genetic_code='1+U') for s in pr.get_fasta(subj_pr, subject)])
    write(f'Protein sequences with Selenocysteine (U)')
    # else:
    #     write(f'Checking protein sequences')
    #     # translates the CDS sequences into protein (no 'U's)
    #     prot_query = pd.Series([translate(s) for s in seq_query])
    #     prot_query.index = query_pr.Q_align_prot_seq.index.copy()
    #     if prot_query.equals(query_pr.Q_align_prot_seq) == True:
    #         write(f'Query protein sequences checked')
    #     else:
    #         # '~' symbol inverse the characters of a series, True is 
    #         # transformed in False and the other way around.
    #         # index = those rows where the 'Q_align_prot_seq' is different
    #         # (False) to the prot_query.
    #         index = ~(prot_query.eq(query_pr.Q_align_prot_seq))
    #         print(seq_query[index])
    #         print(prot_query[index])
    #         print(query_pr.Q_align_prot_seq[index])
    #         raise Exception()

    #     prot_subj = pd.Series([translate(s) for s in seq_subj])
    #     prot_subj.index = subj_pr.Subj_align_prot_seq.index.copy()
    #     if prot_subj.equals(subj_pr.Subj_align_prot_seq) == True:
    #         write(f'Subject protein sequences checked')
    #     else:
    #         index = ~(prot_subj.eq(subj_pr.Subj_align_prot_seq))
    #         print(seq_subj[index])
    #         print(prot_subj[index])
    #         print(subj_pr.Subj_align_prot_seq[index])
    #         raise Exception()
    #     write(f'Tblastx and Pyranges give the same protein sequences')

    query_df = query_pr.as_df() # converts into dataframe
    subj_df = subj_pr.as_df()
    # gc.collect()
    
    return query_df, subj_df

def pandas2(subj_df_gff, query_df_gff, subject, query):
    '''
    This function joins the subject and query gff-format dataframes.

    Parameters
    ----------
    subj_df_gff : Dataframe
        Dataframe with the gff columns related to the subject
    query_df_gff : Dataframe
        Dataframe with the gff columns related to the query
    subject : String, easyterm object
        Path to the file of the subject
    query : String, easyterm object
        Path to the file of the query
    output : String, easyterm object
        Path to the output folder

    Returns
    -------
    joined_df : Dataframe
        Dataframe with all the columns of the tblastx hits
    '''

    query_df, subj_df = translate_prot(subj_df_gff, query_df_gff, subject, query, CDS_sequences=True)
    # rename the PyRanges-format to join both DataFrames
    query_df = query_df.rename(columns={'Chromosome': 'Q_ID', 'Start': 'Q_align_s',
                                        'End': 'Q_align_e', 'Strand': 'Q_Strand'})
    # 'extend_orfs.py' is made to ignore these columns, they are empty
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

    for i in query_df.index:
        # now we need to separate the different columns in 'Attribute'
        for x in query_df.loc[i, 'Attribute'].split(';'):
            if x.startswith('Attribute'):
                # strip function removes white spaces or characters from the 
                # beginning or end of a string.
                list_attribute_query.append(x.split('=')[1].strip())
            elif x.startswith('Q_fr'):
                list_fr_query.append(x.split('=')[1].strip())
        for x in subj_df.loc[i, 'Attribute'].split(';'):
            if x.startswith('Attribute'):
                list_attribute_subj.append(x.split('=')[1].strip())
            elif x.startswith('Value'):
                list_value_subj.append(x.split('=')[1].strip())
            elif x.startswith('Subj_fr'):
                list_fr_subj.append(x.split('=')[1].strip())
            else:
                list_evalue_subj.append(x.split('=')[1].strip())
        # gc.collect()

    query_df['Attribute'] = list_attribute_query
    query_df['Q_fr'] = list_fr_query
    subj_df['Attribute'] = list_attribute_subj
    subj_df['Score'] = list_value_subj
    subj_df['Subj_fr'] = list_fr_subj
    subj_df['Evalue'] = list_evalue_subj

    # merge both dataframes according to 'Attribute' column (old 'ID' column)
    joined_df = subj_df.merge(query_df.set_index(['Attribute']), on='Attribute')
    return joined_df

def pairwise(joined_df, matrix):
    '''
    This function runs pairwise global alignment tool to introduce gaps in the 
    protein sequences.

    Parameters
    ----------
    joined_df : Dataframe
        Dataframe with all the remaining tblastx hits
    matrix : Dictionary of tuples
        Dictionary with the BLOSUM62 matrix values

    Returns
    -------
    joined_df : Dataframe
        Dataframe with all the remaining tblastx hits plus gaps in the protein
        sequences.
    '''

    write(f'Pairwise alignment')
    # for i, row in joined_df.iterrows():
    for i in joined_df.index:
        # -7 is the cost to open a gap, -1 is the cost to extend it
        # pairwise2.align.global parameters:
            # d     A dictionary returns the score of any pair of characters
            # s     Same open and extend gap penalties for both sequences
        alignment = pairwise2.align.globalds(joined_df.at[i, 'Q_align_prot_seq'], 
                                             joined_df.at[i, 'Subj_align_prot_seq'], 
                                             matrix, -7, -1,
                                             one_alignment_only = True)
        # only the best scored alignment is selected
        joined_df.at[i, 'Q_align_prot_seq'] = alignment[0][0]
        joined_df.at[i, 'Subj_align_prot_seq'] = alignment[0][1]
        # gc.collect()

    return joined_df

def UGA(query_prot_seq, subj_prot_seq, u_subj, u_query):
    '''
    Function made to find the 'U' responsible for the readthrough

    Parameters
    ----------
    query_prot_seq : String
        Protein sequence corresponding to the query
    subj_prot_seq : String
        Protein sequence corresponding to the subject
    u_subj : Int
        Index of the first 'U' in the subject protein sequence
    u_query : Int
        Index of the first 'U' in the query protein sequence

    Returns
    -------
    u_subj : Int
        Index of the good 'U' in the subject protein sequence
    u_query : Int
        Index of the good 'U' in the query protein sequence
    '''

    list_align_ugas = list()
    center_alignment_subj = len(subj_prot_seq)/2
    center_alignment_q = len(query_prot_seq)/2

    for idx, x in enumerate(subj_prot_seq):
        if x == 'U' and query_prot_seq[idx] == 'U':
            list_align_ugas.append(idx)
    # maybe the 'U's are not aligned
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

        u_subj = sorted_list[0] # aligned 'U' closer to the center
        u_query = sorted_list[0]

    return u_subj, u_query

def list_UGAs(row, u_subj, u_query):
    '''
    This function creates 4 lists to separate the U's in upstream or downstream
    according to their position regarding the selected selenocysteine 
    (u_subj and u_query).

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
    Second filter of the script. During this filter, we will get only those hits 
    with aligned selenocysteines (U) in query and target protein sequences.

    Parameters
    ----------
    best_hits : Dataframe
        Dataframe with all the remaining hits
    dictionary_matrix : Dictionary of dictionaries
        Dictionary with the BLOSUM matrix values
    conservation_up : Int, easyterm object
        opt['cons_up']
    conservation_down : Int, easyterm object
        opt['cons_down']

    Returns
    -------
    selenoproteins : Dataframe
        Dataframe with the selenoproteins candidates
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
    
    best_hits.drop_duplicates(subset='Query_CDS', #???
                              inplace=True, ignore_index=True)
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
            if u_subj == u_query: # filter for aligned 'U's
                if n_stops_subj >= 2 or n_stops_query >= 2:
                    # updates protein and cds sequences according to whether 
                    # there is U's upstream, downstream or both.
                    (list_up_query, list_down_query, 
                     list_up_subj, list_down_subj) = list_UGAs(row, u_subj, 
                                                               u_query)
                    # we put the downstream first to not change the length of
                    # the sequences when list_up != 0.
                    # cuts the sequences (replacing with gaps) from the closest
                    # stop to the selected 'U'.
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
                # counts the nº of gaps in both sequences
                n_gaps_subj = row['Subj_align_prot_seq'][:u_subj].count('-')
                n_gaps_query = row['Q_align_prot_seq'][:u_query].count('-')
                # when using cds sequences we need to subtrack the number of
                # gaps and multiply by 3 (1 aa = 3 nucleotides).
                index_3t_nucl_subj = (u_subj - n_gaps_subj) * 3
                index_3t_nucl_query = (u_query - n_gaps_query) * 3
                # filters only when the selected 'U' = 'TGA'
                if row['Subj_CDS'][index_3t_nucl_subj:index_3t_nucl_subj + 3] == 'TGA' and (
                        row['Query_CDS'][index_3t_nucl_query:index_3t_nucl_query + 3] == 'TGA'):
                    # measures the conservation values of the alignment, according
                    # to the matrix BLOSUM62 values.
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
                        # measures the new score (if the sequences have changed)
                        row['Score'] = score(row['Q_align_prot_seq'], 
                                             row['Subj_align_prot_seq'], 
                                             dictionary_matrix)
                        list_score.append(row['Score'])
                        # calculates the density score (score/length of the seq)
                        density_score = round(
                            row['Score']/len(row['Q_align_prot_seq']), 4)
                        list_density_score.append(density_score)
                        # putting the iloc in this way, we take the entire row
                        selenoproteins = (
                            selenoproteins.append(best_hits.iloc[[i]], 
                                                  ignore_index=True))
    # updates the columns of the dataframe
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
    return selenoproteins

def conservation(row, dictionary_matrix, u_query):
    '''
    This function measures the score value upstream and downstream of the 
    alignment of query and target.

    Parameters
    ----------
    row : Row of a Dataframe
        Tblastx hit with all the columns
    dictionary_matrix : Dictionary of dictionaries
        Dictionary with the BLOSUM62 matrix values
    u_query : Int
        Index of the selected 'U' in the query protein sequence

    Returns
    -------
    conservation_before_tga : Int
        Conservation score upstream of the alignment
    conservation_after_tga : Int
        Conservation score downstream of the alignment
    '''

    conservation_before_tga = 0
    conservation_after_tga = 0

    for i, x in enumerate(row['Q_align_prot_seq']):
        # gaps and stops are not taken into account to measure score
        if x == '*' or x == '-' or row['Subj_align_prot_seq'][i] == (
                '-' or row['Subj_align_prot_seq'][i] == '*'):
            continue
        # score upstream
        if i < u_query:
            conservation_before_tga += dictionary_matrix[(x, row['Subj_align_prot_seq'][i])]
        # score downstream
        elif i > u_query:
            conservation_after_tga += dictionary_matrix[(x, row['Subj_align_prot_seq'][i])]

    return conservation_before_tga, conservation_after_tga

def make_aligned_cds(subj_aligned_cds, query_aligned_cds, subj_aligned_pep, query_aligned_pep):
    '''
    Function to insert gaps into the nucleotidic sequences of both query
    and subject.

    Parameters
    ----------
    subj_aligned_cds : String
        Nucleotide sequence from the subject
    query_aligned_cds : String
        Nucleotide sequence from the query
    subj_aligned_pep : String
        Protein sequence from the subject
    query_aligned_pep : String
        Protein sequence from the query

    Returns
    -------
    subj_aligned_cds : String
        Nucleotide sequence from the subject, with gaps
    query_aligned_cds : String
        Nucleotide sequence from the query, with gaps
    '''

    # subj_aligned_cds = subj_cds
    # query_aligned_cds = query_cds

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
    downstream of the good U.

    Parameters
    ----------
    row : Row of a Dataframe
        Tblastx hit with all the columns

    Returns
    -------
    u_dN_dS : String/Float
        Non-synonymous/synonymous ratio upstream U
    d_dN_dS : String/Float
        Non-synonymous/synonymous ratio downstream U
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
    # calculates the index of the good U, in case there are two or more 
    # aligned U's.
    u_subj, u_query = UGA(row['Subj_align_prot_seq'], 
                          row['Q_align_prot_seq'],
                          u_subj, u_query)
    # calculates the CDS mutations both up and downstream
    changes_dN_dS_u = count_coding_changes(query_cds[:u_query * 3], 
                                           subj_cds[:u_subj * 3])
    changes_dN_dS_d = count_coding_changes(query_cds[(u_query + 1) * 3:], 
                                           subj_cds[(u_subj + 1) * 3:])
    # calculates the dN/dS ratio both up and downstream
    tupla_dN_dS_u = count_coding_sites(query_cds[:u_query * 3])
    tupla_dN_dS_d = count_coding_sites(query_cds[(u_query + 1) * 3:])
    # c
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
    This function runs blastp, it is used to get the name of the candidates 
    that have passed all the filters.

    Parameters
    ----------
    selenocandidates : Dataframe
        Dataframe with the selenoprotein candidates
    db : Uniprot database, easyterm object
        opt['uniprot']
    n_cpu : Easyterm object
        opt['c']
    blastp_outfile : String
        Path for the outfile of the blastp
    fasta_query_outfile : String
        Path for the fasta file

    Raises
    ------
    Exception
        We raise an exception to stop the program in case returncode returns 
        different than zero, indicating that subprocess.run hasn't run 
        successfully.
        We also print the stdout and stderr to know more about the problem

    Returns
    -------
    blastp_outfile : String
        Path for the outfile of the blastp
    '''

    write(f'Running blastp, results at {blastp_outfile}')

    fasta_outfile = ''
    # creates a fasta file with the query protein sequences
    with open(fasta_query_outfile, 'w') as fw:
        for i, row in selenocandidates.iterrows():
            # blasp alignments must be done with the sequences without the gaps
            q_no_hypens = row['Q_align_prot_seq'].replace('-', '')
            fw.write(f">{row['Q_ID']}\n")
            fw.write(f"{q_no_hypens}\n")
    # command to run blastp with the query, the database (UniRef50) and an outfile
    # -max_hsps option is to select a max nº of hits per query
    format_6_cmd = (
        'blastp -task blastp-fast -query ' + fasta_query_outfile + 
        ' -db ' + db + ' -out ' + blastp_outfile + 
        ' -num_threads ' + str(n_cpu) + ' -max_hsps ' + str(10) + 
        " -outfmt '6 sacc stitle sstart send sframe qacc " + 
        "qstart qend qframe bitscore evalue sseq qseq'")
    # splits the string into a shell-like syntax
    format_6_cmd_list = shlex.split(format_6_cmd)
    # subprocess module allows you to spawn new processes, connect to their 
    # input/output/error pies, and obtain their return codes.
    x = subprocess.run(format_6_cmd_list, capture_output=True) 
    if x.returncode != 0:
        print(x.stderr, x.stdout)
        raise Exception()

    os.remove(fasta_query_outfile)
    # creates a table DataFrame and names the columns
    uniprot_IDs = pd.read_table(blastp_outfile, 
                                names=['Subj_ID', 'Uniprot_Title', 
                                       'Subj_align_s', 'Subj_align_e', 
                                       'Subj_fr', 'Q_ID', 'Q_align_s',
                                       'Q_align_e', 'Q_fr', 'Score', 
                                       'Evalue', 'Subj_align_prot_seq', 
                                       'Q_align_prot_seq'])
    # c
    uniprot_IDs = uniprot_IDs.groupby(
        'Q_ID', as_index = False).agg({'Uniprot_Title': join_titles})
    selenoproteins = selenocandidates.join(
        uniprot_IDs.set_index('Q_ID'), on='Q_ID')

    return selenoproteins, uniprot_IDs

def join_titles(uniprot_list):
    '''
    Generates the Uniprot titles
    
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
    Third and last filter of the script
    
    This function is used to measure the evolutionary conservation between 
    both query and subject protein sequences.
    
    Parameters
    ----------
    selenoproteins_candidates : Dataframe
        Dataframe with all the selenoprotein candidates
    dictionary_sel : List of dictionaries
        List of dictionaries with the values of the BLOSUM Matrix
        
    Returns
    -------
    outfile : String
        Outfile for the comparison file (blast default format)
    selenoproteins : Dataframe
        Final dataframe with the selenoprotein candidates
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
            elif dictionary_sel[(x, row['Subj_align_prot_seq'][index])] > 0:
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
        Dataframe with all the selenoprotein candidates
    cds_q : Boolean, easyterm object
        opt['cds_q']
    path_cds_q : String
        Path to save the cds sequence from the query if opt['cds_q'] == True
    cds_t : Boolean, easyterm object
        opt['cds_t']
    path_cds_t : String
        Path to save the cds sequence from the target if opt['cds_t'] == True
    pep_q : Boolean, easyterm object
        opt['pep_q']
    path_pep_q : String
        Path to save the protein sequence from the query if opt['pep_q'] == True
    pep_t : Boolean, easyterm object
        opt['pep_t']
    path_pep_t : String
        Path to save the protein sequence from the target if opt['pep_t'] == True
    gff_q : Boolean, easyterm object
        opt['dff_q']
    gff_t : Boolean, easyterm object
        opt['dff_t']
    path_query_gff : String
        Path to save the gff file from the query if opt['dff_q'] == True
    path_subj_gff : String
        Path to save the gff file from the target if opt['dff_t'] == True
    dotplot : Boolean, easyterm object
        opt['dotplot']
    path_dotplot : String
        Path to save the dotplot if opt['dotplot'] == True
        
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

def get_proc_status(keys = None):
    '''
    

    Parameters
    ----------
    keys : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    TYPE
        DESCRIPTION.
    '''
    
    with open('/proc/' + str(os.getpid()) + '/status') as f:
        data = dict(map(str.strip, line.split(':', 1)) for line in f)

    return tuple(data[k] for k in keys) if keys else data

def main():

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
    -cds_q : control the creation of the nucleotide sequences of the 
             candidates in the query.
    -cds_t : control the creation of the nucleotide sequences of the 
             candidates in the subject.
    -dotplot : control the creation of one or more dotplots with the 
               candidates sequences.
    -pep_q : control the creation of the protein sequences of the 
             candidates in the query.
    -pep_t : control the creation of the nucleotide sequences of the 
             candidates in the subject.
    -dff_q : control the creation of a .dff file with the annotations of the 
             candidates in the query.
    -dff_t : control the creation of a .dff file with the annotations of the 
             candidates in the subject.
    -UniRef : uniprot blast database
    -cons_up : conservation minimum value upstream
    -cons_down : conservation minimum value downstream
    -n_section : control the rerun of the sections of the scripts

    ### Options:
    -print_opt: print currently active options
    -h | --help: print this help and exit"""

    def_opt= {
        'q': 'query',
        't': 'subject',
        'o': '/users-d3/EGB-invitado4/seleno_prediction/outputs/',
        'c': 4,
        'd': False,
        'f': False,
        'n_chunks': 10,
        'cds_q': False,
        'cds_t': False,
        'dotplot': False,
        'pep_q': False,
        'pep_t': False,
        'gff_q': False,
        'gff_t': False,
        'uniprot': 
            '/users-d3/EGB-invitado4/seleno_prediction/data/uniref50.fasta',
        'cons_up': 50,
        'cons_down': 150,
        'n_section': 1000}

    opt = command_line_options(def_opt, help_msg)

    check_file_presence(opt['q'], 'query')
    # return the base name of pathname path
    # query = path.basename(opt['q']).split('.')[0]
    # subject = path.basename(opt['t']).split('.')[0]

    if not os.path.exists(opt['o']):
        os.makedirs(opt['o'])

    # Main paths
    tblastx_outfile = opt['o'] + 'tblastx.tsv'
    tblastx_df_path = opt['o'] + 'tblastx_df_prechunking.tsv'
    postchunking_file = opt['o'] + 'tblastx_df_postchunking.tsv'
    path_fragments = opt['o'] + 'all_orfs.tsv'
    selected_IDs = opt['o'] + 'nov_orfs.tsv'
    path_table_df = opt['o'] + 'ext_orfs.txt'
    path_pairwise = opt['o'] + 'aln_orfs.tsv'
    path_selenocandidates = opt['o'] + 'candidates.tsv'
    out_blastp = opt['o'] + 'candidates_blastp.tsv'
    comparison_outfile = opt['o'] + 'candidates_pretty.txt'
    # Temporal paths
    query_gff = opt['o'] + 'table_query.gff'
    subject_gff = opt['o'] + 'table_subj.gff'
    out_subj_gff = opt['o'] + 'out_subj.gff'
    out_query_gff = opt['o'] + 'out_query.gff'
    path_fasta_outfile = opt['o'] + 'fasta_seq.txt'
    # Output paths
    path_cds_q = opt['o'] + 'candidates_query.cds.fa'
    path_cds_t = opt['o'] + 'candidates_target.cds.fa'
    path_pep_q = opt['o'] + 'candidates_query.pep.fa'
    path_pep_t = opt['o'] + 'candidates_target.pep.fa'
    path_gff_q = opt['o'] + 'candidates_query.gff'
    path_gff_t = opt['o'] + 'candidates_target.gff'
    path_dot_plot = opt['o'] + 'candidates_dotplot.png'
    # fragments_dot_plot = opt['o'] + 'fragments_dotplot.png'
    # overlapping_dot_plot = opt['o'] + 'overlapping_dotplot.png'

    write(f'TwinStop {__version__}')
    write(f'{current_time}')

    write(f'\n### PHASE 1: TBLASTX')

    # tracemalloc.start()
    run_tblastx(opt['q'], opt['t'], opt['o'], opt['c'], opt['d'], opt['f'], tblastx_outfile)
    if not os.path.exists(postchunking_file) or opt['n_section'] < 2:
        columns = pandas1(tblastx_outfile, tblastx_df_path)
        chunking(tblastx_df_path, opt['n_chunks'], lambda x: conv_pyran(x, opt['t'], opt['q']),
                 columns, postchunking_file)
        df = pd.read_csv(postchunking_file, sep='\t', header=0)
    else:
        write(f'Reading {postchunking_file}')
        df = pd.read_csv(postchunking_file, sep='\t', header=0)
        if len(df) == 0:
            write(f'Empty file {df}')

    dictionary_matrix = dictionary_seleno()
    now_1 = datetime.now()
    time_usage = now_1 - now
    write(f'Time usage: {time_usage}')
    peak_memory, current_memory = get_proc_status(('VmHWM', 'VmRSS'))
    # current_memory, peak_memory = tracemalloc.get_traced_memory()
    write(f'Memory peak: {peak_memory}')
    write(f'Current memory use: {current_memory}')
    # tracemalloc.clear_traces()
    # tracemalloc.reset_peak()

    write(f'\n### PHASE 2: FRAGMENTATION')

    if not os.path.exists(path_fragments) or opt['n_section'] < 3:
        fragments = max_score_block(df, dictionary_matrix)
        fragments.to_csv(sep='\t', path_or_buf = path_fragments)
    else:
        write(f'Reading {path_fragments}')
        fragments = pd.read_csv(path_fragments, sep='\t', header=0, 
                                index_col=0)
        if len(fragments) == 0:
            write(f'Empty file {fragments}')

    del df
    now_2 = datetime.now()
    time_usage = now_2 - now_1
    write(f'Time usage: {time_usage}')
    peak_memory, current_memory = get_proc_status(('VmHWM', 'VmRSS'))
    # current_memory, peak_memory = tracemalloc.get_traced_memory()
    write(f'Memory peak: {peak_memory}')
    write(f'Current memory use: {current_memory}')
    # tracemalloc.clear_traces()
    # tracemalloc.reset_peak()

    # if opt['dotplot'] == True:
    #     dot_plot(fragments, fragments_dot_plot)

    write(f'\n### PHASE 3: OVERLAP FILTER')

    if not os.path.exists(selected_IDs) or opt['n_section'] < 4:
        overlapping(fragments, selected_IDs)
        # clusters.to_csv(sep='\t', path_or_buf = selected_IDs)
        clusters = pd.read_csv(selected_IDs, sep='\t', header=0)
    else:
        write(f'Reading {selected_IDs}')
        clusters = pd.read_csv(selected_IDs, sep='\t', header=0)
        if len(clusters) == 0:
            write(f'Empty file {clusters}')

    del fragments
    # if opt['dotplot'] == True:
    #     dot_plot(clusters, overlapping_dot_plot)

    now_3 = datetime.now()
    time_usage = now_3 - now_2
    write(f'Time usage: {time_usage}')
    peak_memory, current_memory = get_proc_status(('VmHWM', 'VmRSS'))
    # current_memory, peak_memory = tracemalloc.get_traced_memory()
    write(f'Memory peak: {peak_memory}')
    write(f'Current memory use: {current_memory}')
    # tracemalloc.clear_traces()
    # tracemalloc.reset_peak()

    write(f'\n### PHASE 4: EXTEND ORFS')

    if not os.path.exists(path_table_df) or opt['n_section'] < 5:
        table_df = extend_orfs(clusters, subject_gff, query_gff, opt['t'], 
                               opt['q'], out_subj_gff, out_query_gff)
        table_df.to_csv(sep='\t', path_or_buf=path_table_df)
    else:
        write(f'Reading {path_table_df}')
        table_df = pd.read_csv(path_table_df, sep='\t', header=0, index_col=0)
        if len(table_df) == 0:
            write(f'Empty file {table_df}')

    del clusters
    now_4 = datetime.now()
    time_usage = now_4 - now_3
    write(f'Time usage: {time_usage}')
    peak_memory, current_memory = get_proc_status(('VmHWM', 'VmRSS'))
    # current_memory, peak_memory = tracemalloc.get_traced_memory()
    write(f'Memory peak: {peak_memory}')
    write(f'Current memory use: {current_memory}')
    # tracemalloc.clear_traces()
    # tracemalloc.reset_peak()

    write(f'\n### PHASE 5: ALIGNMENTS')

    if not os.path.exists(path_pairwise) or opt['n_section'] < 6:
        joined_df = pairwise(table_df, dictionary_matrix)
        joined_df.to_csv(sep='\t', path_or_buf = path_pairwise)
    else:
        write(f'Reading {path_pairwise}')
        joined_df = pd.read_csv(path_pairwise, sep='\t', header=0, index_col=0)
        if len(joined_df) == 0:
            write(f'Empty file {joined_df}')

    del table_df
    now_5 = datetime.now()
    time_usage = now_5 - now_4
    write(f'Time usage: {time_usage}')
    peak_memory, current_memory = get_proc_status(('VmHWM', 'VmRSS'))
    # current_memory, peak_memory = tracemalloc.get_traced_memory()
    write(f'Memory peak: {peak_memory}')
    write(f'Current memory use: {current_memory}')
    # tracemalloc.clear_traces()
    # tracemalloc.reset_peak()

    write(f'\n### PHASE 6: FILTER')

    if not os.path.exists(path_selenocandidates) or opt['n_section'] < 7:
        candidates = UGA_alignments(joined_df, dictionary_matrix, 
                                    opt['cons_up'], opt['cons_down'])
        candidates.to_csv(sep='\t', path_or_buf = path_selenocandidates)
    else:
        write(f'Reading {path_selenocandidates}')
        candidates = pd.read_csv(path_selenocandidates, 
                                 sep='\t', header=0, index_col=0)
        if len(candidates) == 0:
            write(f'Empty file {selenocandidates}')

    del joined_df
    now_6 = datetime.now()
    time_usage = now_6 - now_5
    write(f'Time usage: {time_usage}')
    peak_memory, current_memory = get_proc_status(('VmHWM', 'VmRSS'))
    # current_memory, peak_memory = tracemalloc.get_traced_memory()
    write(f'Memory peak: {peak_memory}')
    write(f'Current memory use: {current_memory}')
    # tracemalloc.clear_traces()
    # tracemalloc.reset_peak()

    write(f'\n### PHASE 7: BLASTP FOR TITLES')

    if not 'Uniprot_Title' in candidates.columns:
        candidates, candidates_blastp = alignments_blastp(
            candidates, opt['uniprot'], opt['c'], 
            out_blastp, path_fasta_outfile)

    now_7 = datetime.now()
    time_usage = now_7 - now_6
    write(f'Time usage: {time_usage}')
    peak_memory, current_memory = get_proc_status(('VmHWM', 'VmRSS'))
    # current_memory, peak_memory = tracemalloc.get_traced_memory()
    write(f'Memory peak: {peak_memory}')
    write(f'Current memory use: {current_memory}')
    # tracemalloc.clear_traces()
    # tracemalloc.reset_peak()

    write(f'\n### PHASE 8: OUTPUTS')

    candidates_pretty, candidates = evo_conservation(candidates, 
                                                     dictionary_matrix)
    candidates.to_csv(sep='\t', path_or_buf = path_selenocandidates)

    with open(comparison_outfile, 'w') as fw:
        fw.write(candidates_pretty)

    outputs(
        candidates, opt['cds_q'], path_cds_q, opt['cds_t'], path_cds_t,
        opt['pep_q'], path_pep_q, opt['pep_t'], path_pep_t, opt['gff_q'], 
        opt['gff_t'], path_gff_q, path_gff_t, opt['dotplot'], path_dot_plot)

    now_8 = datetime.now()
    time_usage = now_8 - now_7
    write(f'Time usage: {time_usage}')
    peak_memory, current_memory = get_proc_status(('VmHWM', 'VmRSS'))
    # current_memory, peak_memory = tracemalloc.get_traced_memory()
    write(f'Memory peak: {peak_memory}')
    write(f'Current memory use: {current_memory}')
    # tracemalloc.stop()

if __name__ == '__main__':
    main()
