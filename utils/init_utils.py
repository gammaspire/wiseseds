import os
import sys
import re
import numpy as np
from astropy.table import Table
from conversion_utils import clip_negative_outliers, apply_error_floor


###################################################################
# Get list of non-flux Bayes model parameters for PDF diagnostics #
###################################################################
def get_bayes_list(results):
    
    #grab the column names from results.fits
    header_list = results.colnames
    
    #all of the possible wavelength bands
    bands = ['FUV','NUV','WISE1','WISE2','WISE3','WISE4',
         'BASS-g','decamDR1-g','BASS-r','decamDR1-r','decamDR1-z',
            'PACS-blue','PACS-green','PACS-red']
    
    #this one is fun!
    #create a list of all bayes parameters which are NOT fluxes
    #conditions:
        #'bayes' must be in the header label
        #'_err' must not be in the header label
        #the header label cannot be a wavelength band
    
    bayes_list = [x for x in header_list if ('bayes' in x) & ('_err' not in x) & (x.replace('bayes.','') not in bands)]
    
    return bayes_list


################################################################
# Move or delete PDF fits files once used for diagnostic plots #
################################################################
def handle_pdf_fits(destination, galaxy_id, delete_fits=False):
    if delete_fits:
        os.system(f'rm {destination}out/{galaxy_id}*.fits')
    else:
        os.system(f'mv {destination}out/{galaxy_id}*.fits {destination}out/PDF_fits')


#####################################################
# Names of each band in the cigale filter textfile! #
#####################################################
def define_flux_dict(n_or_s):
    
    if n_or_s not in ['n','s']:
        print('Please put "n" for north and "s" for south.')
        sys.exit()
    
    if n_or_s=='n':
        flux_dict_north = {'FUV':'FUV', 'NUV':'NUV', 'G':'BASS-g', 'R':'BASS-r',
                     'W1':'WISE1', 'W2':'WISE2', 'W3':'WISE3', 'W4':'WISE4'}
        return flux_dict_north
    
    flux_dict_south = {'FUV':'FUV', 'NUV':'NUV', 'G':'decamDR1-g', 'R':'decamDR1-r',
                     'Z':'decamDR1-z', 'W1':'WISE1', 'W2':'WISE2', 'W3':'WISE3', 'W4':'WISE4'}
    return flux_dict_south


##############
# Trim Table #
##############
def trim_tables(IDs, redshifts, flux_tab, ext_tab):
    '''
    trim flags according to redshift values (must be positive) and whether the galaxies contain photometry data
    '''
    
    #convert to numpy arrays...
    IDs = np.array(IDs)
    redshifts = np.array(redshifts)
    
    all_flags = (redshifts>0.) # & (flux_tab['photFlag'])
    
    return IDs[all_flags], redshifts[all_flags], flux_tab[all_flags], ext_tab[all_flags]


####################################################
# Check whether dir_path exists...create it if not #
####################################################
def check_dir(dir_path):
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
        print(f'Created {dir_path}')


############################################
# Helper functions to generate input files #
############################################

filter_names_all = ['FUV','NUV','G','R','Z','W1','W2','W3','W4']

def create_fauxtab(params_class, flux_tab, ext_tab, IDs, redshifts):
    
    #isolate needed flux_tab fluxes; convert from nanomaggies to mJy
    #order: FUV, NUV, g, r, (z,) W1, W2, W3, W4
    try:
        flag_n = flux_tab['DEC_MOMENT']>32   #isolates north galaxies
        flag_s = flux_tab['DEC_MOMENT']<32   #isolates south galaxies
    except:
        flag_n = flux_tab['DEC']>32   #isolates north galaxies
        flag_s = flux_tab['DEC']<32   #isolates south galaxies
    
    N=len(flux_tab) #all galaxies...north and south.
    dtype=[('OBJID','str'),('redshift','f4'),('FUV','f4'),('FUV_err','f4'),
           ('NUV','f4'),('NUV_err','f4'),('G','f4'),('G_err','f4'),
          ('R','f4'),('R_err','f4'),('Z','f4'),('Z_err','f4'),
          ('W1','f4'),('W1_err','f4'),('W2','f4'),('W2_err','f4'),
          ('W3','f4'),('W3_err','f4'),('W4','f4'),('W4_err','f4')]
    
    faux_tab = Table(data=np.zeros(N,dtype=dtype))
    faux_tab['OBJID']=IDs
    faux_tab['redshift']=redshifts
        
    #define conversion factor for flux
    conversion_factor = 1.
    #if True, convert fluxes from nanomaggies to mJy
    if params_class.convert_flux:
        conversion_factor = 3.631e-3
        
    #for all filters in filter_names_all...populate the respective data column
    for i in filter_names_all:
        
        fluxes = flux_tab[params_class.flux_id_col + i] * conversion_factor
        flux_errs = flux_tab[params_class.flux_id_col_err + i]   #do not apply conversion factor just yet
        
        ###################
        # Clean the data! #
        ###################
        
        #first create flags to identify every row with no photometry
        no_flux_flag = (fluxes==0.) & (flux_errs==0.)

        #any row with no fluxes will be assigned an np.nan
        fluxes[no_flux_flag] = np.nan
        flux_errs[no_flux_flag] = np.nan
        
        #for galaxies WITH photometry!
        
        #if need to convert invariance to an error...do so
        if params_class.ivar_to_err:
            flux_errs[~no_flux_flag] = np.sqrt(1/flux_errs[~no_flux_flag]) * conversion_factor
            
            #if neither condition is met, convert error as normal
        else:
            flux_errs[~no_flux_flag] = flux_errs[~no_flux_flag] * conversion_factor

        ###################
        # EXTINCTION CORR #
        ###################
        
        #check for flag indicating conversion from transmission to extinction (in magnitudes) is needed
        if params_class.transmission_to_extinction:
            ext_values = -2.5 * np.log10(ext_tab[params_class.extinction_col+i])
        else:
            ext_values = ext_tab[params_class.extinction_col+i]
        
        #Milky Way (MW) extinction corrections (SFD) for each band, given in magnitudes.
        ext_corrections = 10.**(ext_values/2.5)   #converting to linear scale factors
                
        #now apply SFD extinction correction (per Legacy Survey)
        flux_errs[~no_flux_flag] *= ext_corrections[~no_flux_flag]
        fluxes[~no_flux_flag] *= ext_corrections[~no_flux_flag]
        
        #################
        # ERROR FLOORS #
        #################
        
        #If the relative error dF/F < 0.10, then let dF = 0.10*F
        #the idea is that our MINIMUM error floor for fluxes will be set as 10% of the flux value
        #for grz and 15% for W1-4 & NUV+FUV. 
        #CURRENTLY --> using 10% for grz; 13% for FUV, NUV, W1-4.
        flux_errs[~no_flux_flag] = apply_error_floor(fluxes[~no_flux_flag], flux_errs[~no_flux_flag], band=i)

        #####################
        # REMOVING PROBLEMS #
        #####################

        #...another conditional statement. if zero is within the 4-sigma confidence interval of the flux value, keep the negative value. if the flux is OUTSIDE of this limit, then set to NaN. 
        fluxes, flux_errs = clip_negative_outliers(fluxes, flux_errs)
          
        ######################################
        # ADDING FLUX ROWS TO THE FAUX TABLE #
        ######################################

        #appending arrays to faux table
        faux_tab[i] = fluxes
        faux_tab[f'{i}_err'] = flux_errs
        
    faux_tab.add_columns([flag_n,flag_s],names=['flag_north','flag_south'])

    return faux_tab


#generalizing the writing of rows for north and south galaxies...
def write_region(file, table, flux_dict, filter_labels_all, filter_comp_names, region):
    for n in table[table[f'flag_{region}']]:

        #the first two will be ID [0] and redshift [1], by design.
        s_gal = f'{n[0]} {round(n[1],4) } '

        for i in range(len(filter_labels_all)):
            if (filter_comp_names[i]!='Z'):
                if filter_labels_all[i] == flux_dict[filter_comp_names[i]]:
                    flux_val = n[filter_comp_names[i]]
                    flux_err = n[f'{filter_comp_names[i]}_err']
                    #for the flux values and errors, round to 4 decimal places
                    s_gal += f"{flux_val:.4f} {flux_err:.4f} "

                else:
                    s_gal += 'nan nan '
            else:
                s_gal += 'nan nan '

        s_gal += '\n'
        file.write(s_gal)
        
    print(f"{region} galaxies finished", len(table[table[f'flag_{region}']]))

def create_flux_table(params_class, trim=True):
        
    #define flux table, extinction table
    ext_tab = params_class.ext_tab
    flux_tab = params_class.flux_tab
    
    IDs = params_class.IDs
    redshifts = params_class.redshifts
    
    #re-define variables with trimmed data
    if trim:
        IDs, redshifts, flux_tab, ext_tab = trim_tables(IDs, redshifts, flux_tab, ext_tab)
    
    #contains FUV, NUV, G, R, Z, W1, W2, W3, W4, north flag, south flag for all galaxies
    faux_table = create_fauxtab(params_class, flux_tab=flux_tab, ext_tab=ext_tab, IDs=IDs, redshifts=redshifts)
    
    #generate flux dictionaries to map the params.txt flux bands to their CIGALE labels
    flux_dict_north = define_flux_dict('n')
    flux_dict_south = define_flux_dict('s')
    
    bands_north = list(flux_dict_north.keys())
    bands_south = list(flux_dict_south.keys())
    
    #write files...
    check_dir(params_class.dir_path)
        
    with open(params_class.dir_path+'/galaxy_data.txt', 'w') as file:
        
        #create file header!
        s = '# id redshift '
        
        #if the band appears in both hemispheres (like G, R), write only once.
        #otherwise, write labels only if they exist in north/south dict.
        for flux in filter_names_all:
            
            #len(flux)>=2 isolates NUV, FUV, W1-4 bands. excludes gr(z) bands
            if (flux in bands_north) & (flux in bands_south) & (len(flux)>=2):
                s = s + f'{flux_dict_north[flux]} {flux_dict_north[flux]}_err ' #same for N & S
            
            #for G and R, add twice (as we will need one label each for BASS-g north and DECam-g south
            #add labels if needed, else ''.
            #adding them consecutively ensures the file will read something like BASS-g BASS-g_err DECam-g DECam-g_err ...
            #I guess '' is a failsafe in case the band is not in north or south...effectively skips the header label
            #maybe I set that up once upon a time to accommodate Z-band? I don't know.
            else:
                s = s + f'{flux_dict_north[flux]} {flux_dict_north[flux]}_err ' if flux in bands_north else s + ''
                s = s + f'{flux_dict_south[flux]} {flux_dict_south[flux]}_err ' if flux in bands_south else s + ''
            
        s = s + ' \n'
        file.write(s)
        
        #storing header informtaion in a list so I can interate over it
        #however, I am only keeping the CIGALE FILTER LABELS. no #, no id, no redshift
        filter_labels_all = [x for x in s.split() if ('_err' not in x) & (x!='#') & (x!='id') & (x!='redshift')]
        
        #you may wonder why G and R are included TWICE!
            #There are two sources of G and R depending on whether galaxy is in northern
            #or southern hemisphere. The first G (from flux_dict_north) becomes 'BASS-g,' the second 'decamDR1-g'
            #this is how I set up the header file earlier!
        filter_comp_names = ['FUV','NUV','G','G','R','R','Z','W1','W2','W3','W4']
        print(filter_labels_all)
        
        #for every "good" galaxy in flux_tab, add a row to the text file with relevant information
        
        ####################
        ###NORTH GALAXIES###
        ####################
        write_region(file, faux_table, flux_dict_north, filter_labels_all, filter_comp_names, 'north')
        
        ####################
        ###SOUTH GALAXIES###
        ####################
        write_region(file, faux_table, flux_dict_south, filter_labels_all, filter_comp_names, 'south')

        
def create_ini_files(params_class): #dir_path, sfh_module, dust_module, ncores):
    
    check_dir(params_class.dir_path)
    
    #create pcigale.ini files
    with open(params_class.dir_path+'/pcigale.ini', 'w') as file:
        file.write('data_file = galaxy_data.txt \n')
        file.write('parameters_file = \n')
        file.write(f'sed_modules = {params_class.sfh_module}, bc03, nebular, dustatt_modified_CF00, {params_class.dust_module}, skirtor2016, redshifting \n')
        file.write('analysis_method = pdf_analysis \n')
        file.write(f'cores = {params_class.ncores} \n')  

    #create pcigale.ini.spec files
    
    with open(params_class.dir_path+'/pcigale.ini.spec', 'w') as file:
        file.write('data_file = string() \n')
        file.write('parameters_file = string() \n')
        file.write('sed_modules = cigale_string_list() \n')
        file.write('analysis_method = string() \n')
        file.write('cores = integer(min=1)')
    
    
##################################################
# Add a SLEW of module parameters to pcigale.ini #
##################################################
def add_params(dir_path,sed_plots=False,lim_flag='noscaling',nblocks=1, create_pdfs=False):
    
    #different modules, different naming schemes...
    #sfhdelayed --> age_main, age_burst
    #sfh2exp --> age, burst_age

    with open(dir_path+'/pcigale.ini','r') as file:
        lines = file.readlines()
        
    #modify lines...enjoy.
    modified_lines = []
    for line in lines:
        
        if re.match(r'^\s*save_best_sed\s*=', line):
            if sed_plots:
                modified_lines.append("  save_best_sed = True \n")
                print('line changed: save_best_sed = True')
            else:
                modified_lines.append(line)
        
        elif re.match(r'^\s*tau_main\s*=', line):
            modified_lines.append('   tau_main = 300, 500, 1000, 3000, 6000, 1e5 \n')
        
        elif re.match(r'^\s*age\s*=', line):
            modified_lines.append('   age = 1e3, 3e3, 5e3, 7e3, 1e4, 13000 \n') 
        
        elif re.match(r'^\s*age_main\s*=', line):
            modified_lines.append('   age_main = 1e3, 3e3, 5e3, 7e3, 1e4, 13000 \n')  
        
        elif re.match(r'^\s*tau_burst\s*=', line):
            modified_lines.append('   tau_burst = 100, 200, 400 \n')
        
        elif re.match(fr'^\s*burst_age\s*=', line):
            modified_lines.append('   burst_age = 20, 80, 200, 400, 800, 1e3 \n')
        
        elif re.match(fr'^\s*age_burst\s*=', line):
            modified_lines.append('   age_burst = 20, 80, 200, 400, 800, 1e3 \n')    
        
        elif re.match(r'^\s*f_burst\s*=', line):
            modified_lines.append('   f_burst = 0, 0.001, 0.005, 0.01, 0.05, 0.1 \n')
        
        elif re.match(r'^\s*imf\s*=', line):
            modified_lines.append('   imf = 1 \n')
        
        elif re.match(r'^\s*metallicity\s*=', line):
            modified_lines.append('   metallicity = 0.004, 0.02, 0.05 \n')
        
        elif re.match(r'^\s*variables\s*=',line):
            #modified_lines.append('  variables = sfh.sfr, stellar.m_star, sfh.burst_age, sfh.age, sfh.f_burst, sfh.tau_burst, sfh.tau_main, attenuation.Av_ISM, dust.alpha, dust.gamma, dust.qpah, dust.umean, dust.umin, dust.mass \n') 
            modified_lines.append('  variables = sfh.sfr, stellar.m_star, stellar.metallicity, sfh.burst_age, sfh.age, sfh.f_burst, sfh.tau_burst, sfh.tau_main, agn.fracAGN, attenuation.Av_ISM, dust.mass \n') 
        
        elif re.match(r'^\s*normalise\s*=',line):
            modified_lines.append('   normalise = True \n')
        
        elif re.match(r'^\s*Av_ISM\s*=',line):
            modified_lines.append('  Av_ISM = 0.0, 0.01, 0.025, 0.03, 0.035, 0.04, 0.05, 0.06, 0.12, 0.15, 1.0, 1.3, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3 \n')
        
        elif re.match(r'^\s*fracAGN\s*=',line):
            modified_lines.append('  fracAGN = 0.0, 0.05, 0.1, 0.5 \n')
        
        elif re.match(r'^\s*umin\s*=',line):
            modified_lines.append('  umin = 1.0, 5.0, 10.0 \n')
        
        elif re.match(r'^\s*alpha\s*=',line):
            modified_lines.append('  alpha = 1.0, 2.0, 2.8 \n')
        
        elif re.match(r'^\s*gamma\s*=',line):
            modified_lines.append('  gamma = 0.02, 0.1 \n')  
        
        elif re.match(r'^\s*blocks\s*=',line):
            modified_lines.append(f'  blocks = {nblocks} \n')  
        
        elif re.match(r'^\s*lim_flag\s*=',line):
            modified_lines.append(f'  lim_flag = {lim_flag} \n')
        
        elif re.match(r'^\s*save_chi2\s*=',line):
            
            if create_pdfs:
                modified_lines.append("  save_chi2 = properties \n")
                print('line changed: save_chi2 = properties')
            else:
                modified_lines.append(line)            
        else:
            modified_lines.append(line)
    
    #write modified lines back into the file...or, rather, recreate the file with ALL lines, modified or otherwise
    with open(dir_path+'/pcigale.ini', 'w') as file:
        file.writelines(modified_lines)