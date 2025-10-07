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

    
def organize_sed_output(dir_path, main_tab, out_dir_name):
    
    os.chdir(os.path.join(dir_path, out_dir_name))
    print(os.path.join(dir_path, out_dir_name))
    
    print(os.listdir())
    
    #create directory for best SED models
    os.makedirs('best_SED_models', exist_ok=True)
    
    #remove old best model fits files
    os.system('rm *best_model*.fits')
    
    #move best model fits files into folder
    os.system('mv *best_model* best_SED_models')
    
    #remove all SFH fits files
    os.system('rm *_SFH*')