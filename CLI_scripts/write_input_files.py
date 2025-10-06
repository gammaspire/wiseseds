# CLI script for generating the input files for CIGALE

import sys

#covering all bases...just in case.
sys.path.insert(0,'utils')
sys.path.insert(0,'../utils')

from param_utils import Params
from init_utils import create_flux_table, create_ini_files, add_params
from cigale_utils import run_genconf

if __name__ == "__main__":
    
    if '-h' in sys.argv or '--help' in sys.argv:
        print("USAGE: %s [-params <param_file>]")
        sys.exit()

    if '-params' in sys.argv:
        p = sys.argv.index('-params')
        param_file = str(sys.argv[p+1])
    else:
        print('-params argument not found. exiting.')
        sys.exit()

    params = Params(param_file)
    
    #load tables
    params.load_tables()
    
    #load IDs and redshifts!
    params.load_columns()
    
    print('Generating flux table and input .ini files for CIGALE...')
    
    create_flux_table(params)
    create_ini_files(params)
    
    #configure input files and generate configuration files
    print('Configuring input text files...')
    run_genconf(params.dir_path)

    #modify pcigale.ini according to our settings
    add_params(params.dir_path, params.sed_plots, params.lim_flag, params.nblocks, 
               create_pdfs=params.create_pdfs)
    
    print('Input files successfully generated!')