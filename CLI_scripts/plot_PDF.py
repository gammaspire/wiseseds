'''
Standalone CLI (command line interface) script for running the probability distribution function diagnostics
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

#covering all bases...just in case.
sys.path.insert(0,'utils')
sys.path.insert(0,'../utils')
from param_utils import Params
from plotting_utils import generate_pdfs, handle_pdf_fits
from init_utils import get_bayes_list
   
if __name__ == "__main__":

    #unpack params.txt file here
    if '-h' in sys.argv or '--help' in sys.argv:
        print("USAGE: %s [-params (name of parameter.txt file, no single or double quotations marks)]")
    
    if '-params' in sys.argv or '--params' in sys.argv:
        p = sys.argv.index('-params')
        param_file = str(sys.argv[p+1])
    else:
        print('-params not found. exiting.')
        sys.exit()
    
    params = Params(param_file)
    params.find_out()     #define params.output_dir_name. needed!
    
    results = Table.read(f'{params.destination}{params.output_dir_name}/results.fits')
    
    #ensure output directory exists...otherwise, create it!
    pdf_dir = os.path.join(params.destination, params.output_dir_name, 'PDF_fits')
    os.makedirs(pdf_dir, exist_ok=True)
    
    #create PDF
    for index in range(len(results)):
        
        #get galaxy ID
        galaxy_id = results['id'][index]
        
        #get list of bayes parameters!
        bayes_list = get_bayes_list(results)
        
        #generate the PDFs...as indicated by the name of the function
        generate_pdfs(results, params.destination, index, bayes_list, params.output_dir_name)
        
        #removes the .fits used for the PDFs (9 per galaxy)
        handle_pdf_fits(params.destination, params.main_tab, galaxy_id, params.id_col, params.delete_pdf_fits,
                       params.output_dir_name)
        
        print('Finished!')