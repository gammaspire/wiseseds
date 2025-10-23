import sys
import os
from astropy.table import Table, vstack
homedir = os.getenv("HOME")

# ensure utils is in path
sys.path.insert(0, 'utils')
from param_utils import Params   #inherit Params class


def run_cigale_all(herschel=False):
        
    #IF HERSCHEL BANDS, then user must manually complete this following step (e.g., generate 
    #their own .txt files...pending some sort of photometry catalog with row-matched Herschel
    #data).
    if not herschel:
        os.system('python CLI_scripts/write_input_files.py -params params.txt')
    
    os.system('python CLI_scripts/run_cigale_cli.py -params params.txt')

    
if __name__ == "__main__":

    #unpack params.txt file here
    if '-h' in sys.argv or '--help' in sys.argv:
        print("USAGE: %s [-params (name of parameter.txt file, no single or double quotations marks)]")
        sys.exit()
    
    if '-params' in sys.argv:
        p = sys.argv.index('-params')
        param_file = str(sys.argv[p+1])
    else:
        print('-params argument not found. exiting.')
        sys.exit()

    params = Params(param_file)

    herschel = False
    if 'PACS' in params.bands_north:
        herschel = True
                
    run_cigale_all(herschel)
    
    #find most recently edited out*/ directory. NEEDED.
    params.find_out()   #returns self.output_dir_name
    
    if params.create_pdfs:
        
        #read results.fits output
        results_path = os.path.join(params.destination, params.output_dir_name, 'results.fits')
        results = Table.read(results_path)
    
        #ensure output directory exists...
        pdf_dir = os.path.join(params.destination, params.output_dir_name, 'PDF_fits')
        os.makedirs(pdf_dir, exist_ok=True)
        
        #loop through each galaxy to generate PDFs
        for index in range(len(results)):
            os.system('python CLI_scripts/plot_PDF.py -params params.txt')
            
    print(f"Results of this run's SEDs+PDFs (if applicable) and results.fits located in {params.destination}{params.output_dir_name}/.\n")
