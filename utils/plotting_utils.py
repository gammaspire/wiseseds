'''
Various functions for Probability Distribution Function diagnostics.
'''

import warnings
warnings.filterwarnings('ignore')

import sys
import os
import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table

from seaborn import pairplot, load_dataset
import pandas as pd


################################################################
# Move or delete PDF fits files once used for diagnostic plots #
################################################################
def organize_pdf_fits(destination, main_tab, id_col, out_dir_name):
    #move to out/ directory
    os.chdir(os.path.join(destination, out_dir_name))
    
    #create directory for PDF fits
    os.makedirs('PDF_fits', exist_ok=True)
    
    len_tab = len(main_tab)

    #4 for VFID (VFS), 5 for OBJID (WISESize)
    if len(str(id_col)) == 4:
        formatted_strings = [f"{num:04}" for num in range(len_tab)]
        ID_prefix='VFID'
    else:
        formatted_strings = [f"{num:05}" for num in range(len_tab)]
        ID_prefix='OBJID'

    #move all galaxy fits files to PDF_fits
    for num in formatted_strings:
        galID = f'{ID_prefix}{num}'
        os.system(f'mv {galID}*fits PDF_fits')


def handle_pdf_fits(destination, main_tab, galaxy_id, id_col, delete_fits, out_dir_name):
    
    if delete_fits:
        os.system(f'rm {destination}{out_dir_name}/{galaxy_id}*.fits')
    else:
        organize_pdf_fits(destination, main_tab, id_col, out_dir_name)

####################################
# Code to generate the corner plot #
####################################
def corner_plot(df, destination, out_dir_name, galaxy_id, color='orange'):
    
    print(f'PDF finished. generating cornerplot for {galaxy_id}...')
    cplot = pairplot(df,kind='kde',corner=True,diag_kws={'color': color},
            plot_kws={'color': f'dark{color}'})
    cplot.savefig(f'{destination}{out_dir_name}/PDF_fits/{galaxy_id}_corner.pdf')
    plt.close()


#####################
# Generate the PDFs #
#####################
def generate_PDF_plot(results, destination, index, bayes_list, out_dir_name):
    
    #create empty table into which I will add all read-in variables and probabilities
    df = pd.DataFrame([])
    
    galaxy_id = results['id'][index]
    
    try:
        bayes_list.remove('bayes.sfh.tau_main')   #remove; I decided to replace this parameter with agn_frac
        bayes_list.remove('bayes.stellar.metallicity')   #remove also...I enjoy 3x3 grid
    except:
        print('bayes.sfh.tau_main and/or bayes.stellar.metallicity not found. not removed from list of bayes parameters.')
    
    fig, ax = plt.subplots(3, 3,figsize=(26,16))
    fig.suptitle(f'{galaxy_id} Probability Distribution Functions',fontsize=30,y=0.92)
    ax = ax.flatten()
    
    for n, item in enumerate(bayes_list):
        
        x_val = results[item][index]
        x_err = results[item+'_err'][index]
        
        ax[n].axvline(x_val,color='red',ls='-.',label=f'Bayes Value = {x_val:.3e}')
        ax[n].axvspan(x_val-x_err, x_val+x_err, alpha=0.1, color='red')   #shaded region

        item = item.replace('bayes.','')
        prob_tab = Table.read(f'{destination}{out_dir_name}/{galaxy_id}_{item}.fits')

        xcol=[c for c in prob_tab.colnames if c != 'probability'][0]
        
        ax[n].plot(prob_tab[xcol],prob_tab['probability'],marker='o',color='tab:blue')
        
        try:
            #add row to the df
            df[f'probability_{item}'] = list(prob_tab['probability'])
        except Exception as e:
            print('ERROR:', e)
            sys.exit()
        
        best_value = results['best.'+item][index]
        ax[n].axvline(results['best.'+item][index],color='black',ls='--',
                      label=f'Best Value = {best_value:.3e}')

        ax[n].set_xlabel(f'log({item})',fontsize=20)
        ax[n].tick_params(axis='both', labelsize=14)

        ax[n].legend(fontsize=15)

    fig.savefig(f'{destination}{out_dir_name}/PDF_fits/{galaxy_id}_PDF.pdf', 
                bbox_inches='tight', pad_inches=0.2, dpi=100)
    plt.close()
    
    return df, galaxy_id
    

###############################################
# Generate the PDFs AND the corner plot .pdfs #
###############################################
def generate_pdfs(results, destination, index, bayes_list, out_dir_name):
    
    #simple function -- two lines. yay.
    #(or four, if you count these two comments. OR five, if you count the def ____: line)
    df, galaxy_id = generate_PDF_plot(results, destination, index, bayes_list, out_dir_name)
    corner_plot(df, destination, out_dir_name, galaxy_id)