# Crystal_water_analysis
Multistep analysis of crystal structures to identify what positions that are well maintained within the body of structures.

## Requirements:
Python Modules: MDAnalysis
At least AmberTools (really any version). Or any other minimization codes can be substituted.  

## Overview: 
Step 0: Perform a structural alignment of the body of PDB structures. (Not encoded here)

Step 1: Analysis of aligned structures to grab cartesian coordinates of crystal waters as well as b-factor values. Bin the coordinates in a 3D histogram, calculating a frequency and average b-factor for each bin. Weights are included in the config file for this step to allow calculation of a weighted average when there are numerous representations of a single structure.

Following steps are performed for each specific structure:

Step 2: Prepare the structure for parameterization. Allows for setting protonation states of residues as well as renaming residues to fit the parameter set. (Not hugely important in the provided example.) 

Step 3: Run the prepared structure through tleap to add hydrogens to all residues and get a prmtop/inpcrd file for minimization calculations. 

Step 4: Run a basic Sander minimization calculation to allow for hydrogens to optimize positions. Do this step while all waters are present because hydrogen-bonding networks may be drastically changed by removing waters prematurely. 

Step 5: For the given structure, grab positions of all waters. Bin these positions in the calculated 3D histograms and test whether to keep the respective water by testing the 3D histogram values (both frequency and avg. b-factor) with user-defined boolean-tests. 

## To run the example: 
bash automate_docking_struct_prep.sh

## Parameters for dx_analysis.py 

## Parameters for prep_crystal.py 

## Parameters for prep_structures.py 


