# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 11:55:44 2022

@author: Toño
"""

from matplotlib import pyplot as plt
import matplotlib

def dot_plot(final_table, outfile):
    '''
    This function is used to generate one or several dot plots from the 
    information in the pyranges table.
    We define the x and the y using 'loc' method which allows access to a 
    group of rows and columns by label(s).
    And 'iloc' method which allows access to a group of rows and columns 
    by integer position(s).
    
    Parameters
    ----------
    table : pyranges table with the hits of tblastx.
    output : string, path to the folder '/outputs', where we save the results.
    q_name : string, name used for the query in the outfile title.
    s_name : string, name used for the subject in the outfile title.
    
    Returns
    -------
    Dot plot/s are saved in the fileout.
    
    '''
    
    print(f'Making the dot plot')
    
    fig, ax = plt.subplots() # creates two plots
    # contrary to usual python slices, both the start and the stop are included
    x_query = final_table.loc[:, ['Q_align_s', 'Q_align_e']]
    # takes all the rows of the 2 given columns
    y_subject = final_table.loc[:, ['Start', 'End']]
    cmap=matplotlib.cm.get_cmap('viridis') # scale of color
    # define what is the minimum and maximum value of the color gradient
    norm=matplotlib.colors.Normalize(final_table.Score.min(), 
                                     final_table.Score.max())
    scamap= matplotlib.cm.ScalarMappable(norm, cmap) # gradient of color
    
    for i in range(len(final_table)):
        #if final_table['Chromosome'][i] == 'TRINITY_DN1738_c0_g1_i1' and final_table['Q_ID'][i] == 'TRINITY_DN583_c0_g1_i1':
            if int(final_table['Q_fr'][i]) <= 0:
                # takes the 'i' row of the first column: 'Q_align_s'
                x2 = x_query.iloc[i,0] 
                x1 = x_query.iloc[i,1]
            else:
                x1 = x_query.iloc[i,0]
                x2 = x_query.iloc[i,1]
                
            if int(final_table['Subj_fr'][i]) <= 0:
                y2 = y_subject.iloc[i,0]
                y1 = y_subject.iloc[i,1]
            else:
                y1 = y_subject.iloc[i,0]
                y2 = y_subject.iloc[i,1]
            # defines the linestyle depending on the frame
            if final_table['Q_fr'][i] == 1:
                line_style = 'solid'
            elif final_table['Q_fr'][i] == 2:
                line_style = 'dotted'
            else:
                line_style = 'dashed'
                
            dx = x2 - x1
            dy = y2 - y1
            # creates an arrow without the head
            plt.arrow(x1, y1, dx, dy, length_includes_head=True, width=0,
                      color=scamap.to_rgba(final_table.Score[i]), 
                      linestyle=line_style, alpha=0.8)
            
            if dx < 0:
                dx = -0.1
            else:
                dx = 0.1 
            if dy < 0:
                dy = -0.1
            else:
                dy = 0.1
                
            plt.arrow(x2 - dx, y2 - dy, dx, dy, length_includes_head=True, 
                      width=0, head_width = 3,
                      color=scamap.to_rgba(
                          final_table.Score[i]), alpha=0.8) # creates only a head
        
    plt.xlabel(f"{final_table['Q_ID'][i]}")
    plt.ylabel(f"{final_table['Chromosome'][i]}")
    fig.colorbar(scamap, ax=ax) # colorbar of the gradient
    plt.savefig(outfile, dpi=1200)
        
        