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
from plotting_utils import generate_pdfs
   
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
    
    results = Table.read(f'{params.destination}out/results.fits')
    
    #ensure output directory exists...otherwise, create it!
    pdf_dir = os.path.join(params.destination, 'out', 'PDF_fits')
    os.makedirs(pdf_dir, exist_ok=True)
    
    #create PDF
    for index in range(len(results)):
        generate_pdfs(results, params.destination, index, params.delete_pdf_fits)
        print('Finished!')