# CLI script to initialize and execute CIGALE!

import sys

#covering all bases...just in case.
sys.path.insert(0,'utils')
sys.path.insert(0,'../utils')

from param_utils import Params
from cigale_utils import run_cigale, run_sed_plots, organize_sed_output

if __name__ == "__main__":

    #help message
    if '-h' in sys.argv or '--help' in sys.argv:
        print("USAGE: %s [-params (name of parameter.txt file, no single or double quotations marks)]")
        sys.exit()

    #check for params file
    if '-params' in sys.argv:
        p = sys.argv.index('-params')
        param_file = str(sys.argv[p+1])
    else:
        print('-params argument not found. exiting.')
        sys.exit()

    #define ALL parameters using Params class (utils/param_utils.py)
    params = Params(param_file)

    #in order to save the probability distribution functions, ncores = nblocks = 1
    #note that ncores = nblocks = 1 is already established in param_utils.py if params.create_pdfs=1.
    #this print statement is just a confirmation message.
    if params.create_pdfs:
        print('Create PDFs set to True! nblocks = ncores = 1.')

    #run CIGALE
    print('Executing CIGALE...')
    run_cigale(params.destination)

    params.find_out()     #determine most recently edited out*/ directory. needed!
    
    #if SED plots requested, generate them and organize output
    if params.sed_plots:
        print('Generating SED plots...')
        run_sed_plots(params.destination)

        print('Organizing output...')
        organize_sed_output(params.destination, params.main_tab, params.output_dir_name)

    print('CIGALE is Fin!')