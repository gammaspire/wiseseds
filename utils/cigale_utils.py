import os

##################################
# Run CIGALE -- helper functions #
##################################

def run_genconf(dir_path):
    os.chdir(dir_path)
    #will generate configuration files pcigale.ini and pcigale.ini.spec
    #can check/edit the various parameters in the file before the next step
    os.system('pcigale genconf')

    
def run_cigale(dir_path):

    os.chdir(dir_path)   #just in case...
    
    #once completed, this line will generate an 'out' directory containing 'results.fits', among other files.
    os.system('pcigale run')

    
def run_sed_plots(dir_path):
    
    os.chdir(dir_path)  #AGAIN, just in case...
    os.system('pcigale-plots sed')
    print(f'Find .pdf plots in {dir_path}/out/')

    
def organize_sed_output(dir_path, main_tab, id_col):
    os.chdir(os.path.join(dir_path, 'out'))
    
    #create directory for best SED models
    os.makedirs('best_SED_models', exist_ok=True)
    
    #remove old best model fits files
    os.system('rm *best_model*.fits')
    
    #move best model fits files into folder
    os.system('mv *best_model* best_SED_models')
    
    #remove all SFH fits files
    os.system('rm *_SFH*')

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