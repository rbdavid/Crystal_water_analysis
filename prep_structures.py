# ----------------------------------------
# USAGE:
# ----------------------------------------
# python3 prep_structures.py prep_structures.config IO.py

# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys
import importlib
import MDAnalysis
from MDAnalysis.analysis.align import *
import numpy as np

# ----------------------------------------
# VARIABLE DECLARATION
# ----------------------------------------

config_file = sys.argv[1]
IO_functions_file = sys.argv[2]

config_parser = importlib.import_module(IO_functions_file.split('.py')[0],package=None).structure_config_parser
summary = importlib.import_module(IO_functions_file.split('.py')[0],package=None).summary

# ----------------------------------------
# FUNCTIONS: 
# ----------------------------------------

def main():
    # ----------------------------------------
    # LOAD IN THE STRUCTURE TO BE ANALYZED 
    # ----------------------------------------
    u = MDAnalysis.Universe(parameters['prmtop'],parameters['nc'])
    u.trajectory[-1]
    u_all = u.select_atoms('all')
    u_align = u.select_atoms(parameters['mobile_landmark'])
    selection = u.select_atoms(parameters['selection'])
    u_all.translate(-u_align.center_of_mass())

    # ----------------------------------------
    # LOAD IN THE ALIGNMENT STRUCTURE AND ALIGN 
    # ----------------------------------------
    target = MDAnalysis.Universe(parameters['alignment_pdb'])
    target_all = target.select_atoms('all')
    target_align = target.select_atoms(parameters['alignment_landmark'])
    print(target_align.center_of_mass())
    target_all.translate(-target_align.center_of_mass())
   
    print(target_align.n_atoms,u_align.n_atoms)
    R, rmsd = rotation_matrix(u_align.positions,target_align.positions)
    u_all.rotate(R)
    u_all.write('temp.pdb')

    # ----------------------------------------
    # LOAD IN THE DX FILES TO BE ANALYZED
    # ----------------------------------------
    dx_data_container = []
    edges_container = []
    for dx_file in parameters['dx_files']:
        dx_data = []
        with open(dx_file,'r') as f:
            nBins = tuple([int(i) for i in f.readline().split()[-3:]])
            origin = [float(i) for i in f.readline().split()[-3:]]
            delta = [float(f.readline().split()[1]),float(f.readline().split()[2]),float(f.readline().split()[3])]
            test_nBins = tuple([int(i) for i in f.readline().split()[-3:]])
            nElements = int(f.readline().split()[-3])
            if nBins != test_nBins or nBins[0]*nBins[1]*nBins[2] != nElements:
                print('Something with %s is messed up. Killing job.'%(dx_file))
                sys.exit()
            
            #real_origin = [origin[0],origin[1],origin[2]]
            # fixing the origin values because of manipulation during creation of the dx file; VMD weirdness while drawing isosurfaces...
            real_origin = [origin[0]-delta[0]/2.,origin[1]-delta[1]/2.,origin[2]-delta[2]/2.]
            
            x_edges = np.array([real_origin[0]+i*delta[0] for i in range(nBins[0]+1)])
            y_edges = np.array([real_origin[1]+i*delta[1] for i in range(nBins[1]+1)])
            z_edges = np.array([real_origin[2]+i*delta[2] for i in range(nBins[2]+1)])

            for line in f:
                try:
                    dx_data.extend([float(i) for i in line.split()])
                except: # would normally error out due to non-numerical values being turned into floats.
                    continue
            dx_data = np.array(dx_data)     # numpy array of length nElements
            dx_data = np.reshape(dx_data,nBins) # numpy 3d array of size nBins

            dx_data_container.append(dx_data)
            edges_container.append([x_edges,y_edges,z_edges])

    # ----------------------------------------
    # COMPARE WATERS IN STRUCTURE TO 3D HISTOGRAMS
    # ----------------------------------------
    resid_set = []
    for i in range(len(dx_data_container)):
        data = dx_data_container[i]
        x_edges = edges_container[i][0]
        y_edges = edges_container[i][1]
        z_edges = edges_container[i][2]
        delta = (x_edges[1] - x_edges[0],y_edges[1] - y_edges[0],z_edges[1] - z_edges[0])
        max_value = np.max(data)
        print(delta,max_value)
        for j in range(selection.n_residues):
            pos = selection.residues[j].atoms[0].position   # water oxygen is first atom in residue
            x_index = int((pos[0] - x_edges[0])/(delta[0]))
            y_index = int((pos[1] - y_edges[0])/(delta[1]))
            z_index = int((pos[2] - z_edges[0])/(delta[2]))
            try:
                bin_value = data[x_index][y_index][z_index]
                if not eval(parameters['dx_booleans'][i]):  # if booleans return False, return True; add residue to be thrown away...
                    resid_set.extend(selection.residues[j].resid)
            except:
                resid_set.append(selection.residues[j].resid)
    
    resid_set = sorted(list(set(resid_set)))
    keep_selection = 'not (' + parameters['selection'] + ' and resid'
    for i in resid_set:
        keep_selection += ' %s'%(i)
    keep_selection += ')'
    print(keep_selection)
    output_system = u.select_atoms(keep_selection)
    output_system.write(parameters['output_file_name'])
    
    # ----------------------------------------
    # OUTPUTTING SUMMARY INFORMATION
    # ----------------------------------------
    if parameters['summary_boolean']:
        summary('summary.txt',sys.argv,parameters)

# ----------------------------------------
# CREATING PARAMETER DICTIONARY
# ----------------------------------------
parameters = {}
config_parser(config_file,parameters)

# ----------------------------------------
# MAIN
# ----------------------------------------
if __name__ == '__main__':
    main()

