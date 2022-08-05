#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 15:04:27 2022

@author: Toño
"""

import sys
sys.path.insert(0, '/home/mmariotti/bioinfo_dev/pyranges/pyranges')
import pyranges as pr, pandas as pd
from easyterm import command_line_options

def main():
    
    help_msg = '''This program allows us to calculate how many selenoproteins 
    were predicted by both Selenoprofiles and Twinstop and how many did not.
    
    ### Input:
    -selenoprofiles_gff : query gff file predicted by Selenoprofiles
    -twinstop_gff : query gff file predicted by Twinstop
    
    ### Options:
    -print_opt: print currently active options
    -h | --help: print this help and exit'''
    
    def_opt= {
        'selenoprofiles_gff': '/home/ssanchez/seleno_prediction/outputs/',
        'twinstop_gff': '/home/ssanchez/seleno_prediction/outputs/'}
    
    opt = command_line_options(def_opt, help_msg)
    
    #################
    ### load selenoprofiles ref 
    d = pd.read_csv(opt['selenoprofiles_gff'], 
                    sep='\t', header=None,
                    names="Chromosome Source Feature Start End Score Strand x Gene_ID".split())
    
    d = d[ d.Feature=='CDS' ].drop('Source Feature Score x'.split(), axis=1)
    d.Start -= 1 
    
    # predictions may have more than one exon -- makes no sense for 
    # transcriptome.
    # Fixing it by taking boundaries
    dr = pr.PyRanges(d).boundaries('Gene_ID') #, agg={'label':'first', 'Uniprot_Title':'first'})
    dr.label = 'selenoprofiles'
    dr.Uniprot_Title='/'
    
    # checking that preditions don't overlap among each other
    assert((dr.cluster().df.groupby('Cluster').Gene_ID.count() == 1).all())
    
    #################
    ### load twinstop
    tr = pr.read_gff3(opt['twinstop_gff']).drop(
        'Source Feature Frame'.split())
    tr.label = 'twinstop'
    
    before_rem_red=len(tr)
    
    # there may be predictions overlapping! Removing redundancy: take best scoring one
    ## Sergio: you should have a column with Score here! 
    # now it's choosing a random prediction among those overlapping since there's no score to help choosing
    
    tr= pr.PyRanges(tr.cluster().df.sort_values('Score', ascending=False).groupby(
        'Cluster').first().reset_index(drop=True)).drop('Score')
    
    print(f'Remove redundancy : from {before_rem_red} predictions to {len(tr)}\n')
    
    print(f'Ready for benchmark\nSelenoprofiles predictions: {len(dr)}\nTwinstop predictions: {len(tr)}\n')
    bothr = pr.PyRanges(pd.concat([dr.df, tr.df])).cluster()
    
    # count genes per cluster
    c = bothr.df.groupby('Cluster').Gene_ID.count()
    
    # get those lines in clusters with >1 gene
    overlap = bothr[bothr.Cluster.isin(c[c>1].index)] 
    
    # checking that every cluster with >1 predictions has exactly 2
    assert (overlap.df.groupby('Cluster').Gene_ID.count() == 2).all()
    
    tp = overlap[overlap.label == 'twinstop']
    fp = tr[ ~tr.Gene_ID.isin(tp.Gene_ID) ]
    fn = dr[ ~dr.Gene_ID.isin(overlap[overlap.label == 'selenoprofiles'].Gene_ID)]
    
    print(f'True Positives: {len(tp)}\nFalse Positives: {len(fp)}\nFalse Negatives: {len(fn)}')
    
if __name__ == '__main__':
    main()
   
