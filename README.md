# WISESEDS -- Wrapper Scripts to Generate Spectral Energy Distribution Models for WISESize Galaxies

See the Wiki for this code [here](https://github.com/gammaspire/wisesize/wiki/CIGALE-Initialization-and-Execution#cigale-setup) for detailed setup instructions!

## Run the main script:

To run (after navigating to this directory locally):
```
conda activate cigale
python run_cigale.py -params params.txt
```
This main script will run CIGALE on both north and south galaxies, with specific parameters according to params.txt. Output (and input) files will default to specific directories as indicated in params.txt.

To edit the parameter ranges written in the pcigale.ini file (which CIGALE interfaces with directly -- the params.txt file simply streamlines the automatic creation of pcigale.ini), then please refer to `$utils/init_utils.py`. Note that the loop is set up such that parameters which none of the modules use will NOT be included in the mature pcigale.ini file...so try not to edit irrelevant parameter ranges.

The current setup of this code facilitates the running of CIGALE for a parent sample of galaxies belonging to the Northern and Southern hemispheres per their (declination) coordinates in RA-DEC space. I have tried to generalize the script as much as possible, but I have tested it on but two sets of galaxy logging conventions (VFS, WISESize) and some of the labeling may still remain too specific to these conventions (e.g., column names). In these cases where naming errors occur, try poking and prodding the `write_input_files.py` script.

### A CAUTIONARY NOTE: If you will be including Herschel bands and this data are NOT row-matched to the main tables, then the user must manually generate their own photometry.txt files.

## /utils
- A slew of general helper functions (basically organized) that the scipts will call for parsing the params.txt file, setting up and running CIGALE, and optionally generating PDF plots.

## /CLI_scripts
- Standalone scripts that can be run individually as a command line (literally, Command Line Interface).
    - write_input_files.py -- will output, in the directory indicated in params.txt, the files needed to initialize CIGALE. These include pcigale.ini, pcigale.ini.spec, and galaxy_data.txt (photometry tables written in a CIGALE-friendly format).
    - run_cigale_cli.py -- the CLI to run CIGALE, assuming `write_input_files.py` has already been executed. If the user has marked SED_plots=1 in params.txt, then running this script will also generate these SED figures. This script will NOT generate PDF figures.
    - plot_PDF.py -- will generate probability distribution function diagnostics. ee the [Wiki](https://github.com/gammaspire/wisesize/wiki/CIGALE-Initialization-and-Execution#cigale-setup) for instructions. If you need a .diff file, please contact me!
     
## /pcigale_ini_examples
- Two examples of how a mature pcigale.ini and pcigale.ini.spec will look.

## /notebooks
- Jupyter notebook(s), largely for my own purposes of evaluating output CIGALE data.