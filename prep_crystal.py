# ----------------------------------------
# USAGE:
# ----------------------------------------
# python3 prep_crystal.py prep_crystal.config IO.py

# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys
import importlib
import MDAnalysis
import numpy as np

# ----------------------------------------
# VARIABLE DECLARATION
# ----------------------------------------

config_file = sys.argv[1]
IO_functions_file = sys.argv[2]

config_parser = importlib.import_module(IO_functions_file.split('.py')[0],package=None).crystal_config_parser

# ----------------------------------------
# CREATING PARAMETER DICTIONARY
# ----------------------------------------
parameters = {}
config_parser(config_file,parameters)

# ----------------------------------------
# MAIN
# ----------------------------------------
u = MDAnalysis.Universe(parameters['pdb'])
print('Read in %s'%(parameters['pdb']))
u_sel = u.select_atoms(parameters['selection'])
print('Created selection "%s" which has %d atoms.'%(parameters['selection'],u_sel.n_atoms))

for res in parameters['change_resnames']:
    temp_res_index = np.where(u_sel.residues.resids == int(res[0]))[0][0]
    resname = u_sel.residues[temp_res_index].resname
    u_sel.residues[temp_res_index].resname = res[1]
    print('Changed resid %s%s to %s%s'%(resname,res[0],res[1],res[0]))

u_sel.write(parameters['output_file_name'])
print('Saved this selection to %s'%(parameters['output_file_name']))

