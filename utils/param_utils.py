#create dictionary with keyword and values from param textfile

from astropy.table import Table, Row
from conversion_utils import get_redshift


#used for reading the params.txt file!
#returns a dictionary object
def read_params(param_file):

    param_dict={}
    with open(param_file) as f:
        for line in f:
            try:
                key = line.split()[0]
                val = line.split()[1]
                param_dict[key] = val
            except:
                continue
                
    return param_dict


#define a class...easier for me to organize parameters!
class Params():
    
    #################################################
    # extract parameters and assign to variables... #
    #################################################
    
    def __init__(self, param_file):
        
        param_dict = read_params(param_file)

        self.path_to_repos = param_dict['path_to_repos']
        self.main_table = param_dict['main_table']
        self.phot_table = param_dict['phot_table']
        self.extinction_table = param_dict['extinction_table']

        self.dir_path = param_dict['destination']
        self.destination = param_dict['destination']   #I use both interchangeably. do not ask me to fish around 
                                                       #for inconsistencies - I will not.

        self.id_col = param_dict['galaxy_ID_col']

        self.bands_north = param_dict['bands_north'].split("-")   #list!
        self.bands_south = param_dict['bands_south'].split("-")   #list!

        self.ncores = param_dict['ncores']
        self.nblocks = param_dict['nblocks']

        #in order to save the probability distribution functions, ncores = nblocks = 1
        self.create_pdfs = bool(int(param_dict['create_pdfs']))
        self.delete_pdf_fits = bool(int(param_dict['delete_PDF_fits']))
        if self.create_pdfs:
            self.ncores = 1
            self.nblocks = 1

        self.lim_flag = param_dict['lim_flag']
        
        self.sfh_module = param_dict['sfh_module']
        self.dust_module = param_dict['dust_module']
        
        self.Vcosmic_column = param_dict['Vcosmic_column']
        self.redshift_column = param_dict['redshift_column']
        
        self.flux_id_col = param_dict['flux_ID_col']
        self.flux_id_col_err = param_dict['flux_ID_col_err']
        
        self.extinction_col = param_dict['extinction_col']
        
        self.convert_flux = bool(int(param_dict['nanomaggies_to_mJy']))
        self.ivar_to_err = bool(int(param_dict['IVAR_to_ERR']))
        self.transmission_to_extinction = bool(int(param_dict['transmission_to_extinction']))
        
        self.sed_plots = bool(int(param_dict['sed_plots']))
        
        self.load_tables()
        
    ##################################################
    # class functions for loading tables and columns #
    ##################################################

    def load_tables(self):
        #suppress warning text...do not want, do not need.
        import warnings
        from astropy.units import UnitsWarning
        warnings.filterwarnings("ignore", category=UnitsWarning)
                
        #load the tables
        self.main_tab = Table.read(self.path_to_repos + self.main_table)
        self.flux_tab = Table.read(self.path_to_repos + self.phot_table)
        self.ext_tab = Table.read(self.path_to_repos + self.extinction_table)
            
    def load_columns(self):
        self.IDs = self.main_tab[self.id_col]

        #if "vf" in the vcosmic_table name, then must be using Virgo catalogs...thus, Vcosmic column is available
        if 'vf' in self.phot_table:
            Vcosmic_array = self.main_tab[self.Vcosmic_column]
            self.redshifts = get_redshift(Vcosmic_array)

        #otherwise, just use the redshift column
        else:
            self.redshifts = self.main_tab[self.redshift_column]
            
        
    #needed to find the most recent output directory to which CIGALE is sending the results of its most recent run.
    def find_out(self):
        import glob
        import os
        import numpy as np

        pattern = os.path.join(self.destination, "*_out")   #all CIGALE directories end with "_out"
        out_dirs = glob.glob(pattern)
        if not out_dirs or len(out_dirs) < 2:      # if none are found...well. too bad.
            #print(f"No out*/ directories found in {self.destination}. Defaulting to out/.")
            self.output_dir_name = 'out'
            return

        #convert YYYYMMDD_HHMMSS part of the directory name to integers and compare!
        #can do this directly with np.sort(), then pull the last directory name in the array
        latest_out = np.sort(out_dirs)[-1]

        #isolate just the directory name. :-)
        self.output_dir_name = os.path.basename(latest_out)
        print(self.output_dir_name)