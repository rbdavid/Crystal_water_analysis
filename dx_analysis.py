
# ----------------------------------------
# USAGE:
# ----------------------------------------
# python3 dx_analysis.py config_file.config IO.py

# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys
import os
import importlib
import numpy as np
import MDAnalysis

# ----------------------------------------
# VARIABLE DECLARATION: 
# ----------------------------------------

config_file = sys.argv[1]
IO_functions_file = sys.argv[2]

config_parser = importlib.import_module(IO_functions_file.split('.py')[0],package=None).dx_config_parser
summary = importlib.import_module(IO_functions_file.split('.py')[0],package=None).summary

# ----------------------------------------
# FUNCTIONS: 
# ----------------------------------------

def main():
    # ----------------------------------------
    # FILE NAMING VARIABLES
    # ----------------------------------------
    summary_file_name = parameters['output_directory'] + 'summary.txt'
    dx_file_name = parameters['output_directory'] + parameters['dx_output_file_name']
    if parameters['other_cv']:
        cv_file_name = parameters['output_directory'] + parameters['cv_file_name'] 
    # ----------------------------------------
    # LOAD IN THE STRUCTURES TO BE ANALYZED 
    # ----------------------------------------
    positions = np.array([[0.,0.,0.]])
    weights = np.array([0.])
    if parameters['other_cv']:
        cvs = np.array([0.])

    for pdb, weight in parameters['pdb_list']:
        u = MDAnalysis.Universe(pdb)
        selection = u.select_atoms(parameters['selection'])
        weights = np.append(weights,[weight for i in range(selection.n_atoms)])
        positions = np.append(positions,selection.positions,axis=0)
        if parameters['other_cv']:
            cvs = np.append(cvs,calc_cv(selection))
    
    # ----------------------------------------
    # ANALYZE POSITION DATA AND CREATE DX FILES
    # ----------------------------------------
    if parameters['other_cv']:
        half_max_counts = create_dx(positions,weights,float(parameters['dx_delta']),dx_file_name,cv_data=cvs,cv_file_name=cv_file_name)
    else:
        half_max_counts = create_dx(positions,weights,float(parameters['dx_delta']),dx_file_name)
    
    # ----------------------------------------
    # OUTPUTTING SUMMARY INFORMATION
    # ----------------------------------------
    if parameters['summary_boolean']:
        summary(summary_file_name,sys.argv,parameters)

# ----------------------------------------
# CREATING PARAMETER DICTIONARY
# ----------------------------------------
parameters = {}
config_parser(config_file,parameters)

# ----------------------------------------
# LOADING IN NECESSARY FUNCTIONS FROM MODULE FILES
# ----------------------------------------

create_dx = importlib.import_module(parameters['visualization_functions_file'].split('.')[0],package=None).create_dx
if parameters['other_cv']:
    if parameters['cv_type'].upper() == "B_FACTOR":
        calc_cv = importlib.import_module(parameters['user_functions_file'].split('.')[0],package=None).bfactors
    else:
        print('User needs to read in an accepted "cv_type" or edit the main code (line 74) to accept their new function.')
        sys.exit()

# ----------------------------------------
# PREP OUTPUT DIRECTORY STEP
# ----------------------------------------
if parameters['output_directory'][-1] != os.sep:
    parameters['output_directory'] += os.sep

if os.path.exists(parameters['output_directory']):
    print('The output directory, ', parameters['output_directory'], 'already exists. Please select a different directory name for output.')
    sys.exit()
else:
    os.mkdir(parameters['output_directory'])

# ----------------------------------------
# MAIN
# ----------------------------------------
if __name__ == '__main__':
    main()

