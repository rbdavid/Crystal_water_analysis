
# PREAMBLE:

import numpy as np

def create_dx(node_traj_data,node_traj_weights,delta,dx_file_name,cv_data=[],cv_file_name=None):
    # binning node_traj_data that is assumed to be cartesian coordinates.
    x_minmax = (np.min(node_traj_data[:,0]),np.max(node_traj_data[:,0]))
    y_minmax = (np.min(node_traj_data[:,1]),np.max(node_traj_data[:,1]))
    z_minmax = (np.min(node_traj_data[:,2]),np.max(node_traj_data[:,2]))
  
    # calc'ing the desired number of bins for each dimension; adding one bin beyond the max and min for visualization purposes. also, avoids issues with binning
    x_bins = int((x_minmax[1] - x_minmax[0] + 2*delta)/delta)
    y_bins = int((y_minmax[1] - y_minmax[0] + 2*delta)/delta)
    z_bins = int((z_minmax[1] - z_minmax[0] + 2*delta)/delta)

    # creating edges; one bin beyond the min and max values as described above
    x_edges = np.linspace(x_minmax[0]-delta,x_minmax[1]+delta,num=x_bins)
    y_edges = np.linspace(y_minmax[0]-delta,y_minmax[1]+delta,num=y_bins)
    z_edges = np.linspace(z_minmax[0]-delta,z_minmax[1]+delta,num=z_bins)

    nBins = x_bins * y_bins * z_bins        # total number of bins
    
    # prepping binning variables
    binned_space = np.zeros((x_bins,y_bins,z_bins),dtype=np.float64)
    weighted_binned_space = np.zeros((x_bins,y_bins,z_bins),dtype=np.float64)
    if cv_data != []:
        binned_cv = np.zeros((x_bins,y_bins,z_bins),dtype=np.float64)
    
    # looping through all traj data and binning in 3D space
    for i in range(len(node_traj_data)):
        pos = node_traj_data[i]
        x_index = int((pos[0] - x_edges[0])/delta)
        y_index = int((pos[1] - y_edges[0])/delta)
        z_index = int((pos[2] - z_edges[0])/delta)
        
        binned_space[x_index,y_index,z_index] += 1
        weighted_binned_space[x_index,y_index,z_index] += node_traj_weights[i]
        
        # if user has read in cv_data, then calculate the sum of cv data
        if cv_data != []:
            binned_cv[x_index,y_index,z_index] += cv_data[i]
    
    # organization steps to prep for dx output
    flattened_space = binned_space.flatten()
    weighted_flattened_space = weighted_binned_space.flatten()
    half_max_counts = np.max(weighted_flattened_space)/2.0   # can be used for vis state settings
   
    # finishing the average by dividing sum of cv data by the counts; will likely have a large amount of divide by zero warnings
    if cv_data != []:
        flattened_cv = binned_cv.flatten()
        flattened_cv /= flattened_space
        # need to remove a ton of nans
        where_are_nans = np.isnan(flattened_cv)
        flattened_cv[where_are_nans] = 0

    # outputting dx files, which have finicky formats. 
    if nBins%3 != 0:
        # deal with space histogram
        print_able_data = weighted_flattened_space[:-(nBins%3)]
        reshaped = np.reshape(print_able_data, (int(len(weighted_flattened_space)/3),3))
        append_able_data = weighted_flattened_space[-(nBins%3):]
        ## origin not shifted by half delta...
        #np.savetxt(dx_file_name,reshaped,comments='',header='object 1 class gridpositions counts %d %d %d\norigin %.18e %.18e %.18e\ndelta %.18e 0 0\ndelta 0 %.18e 0\ndelta 0 0 %.18e\nobject 2 class gridconnections counts %d %d %d\nobject 3 class array type double rank 0 items %d data follows'%(x_bins,y_bins,z_bins,x_edges[0],y_edges[0],z_edges[0],delta,delta,delta,x_bins,y_bins,z_bins,nBins))    #,footer='object "density (all) [A^-3]" class field'
        # origin shifted by half delta...
        np.savetxt(dx_file_name,reshaped,comments='',header='object 1 class gridpositions counts %d %d %d\norigin %.18e %.18e %.18e\ndelta %.18e 0 0\ndelta 0 %.18e 0\ndelta 0 0 %.18e\nobject 2 class gridconnections counts %d %d %d\nobject 3 class array type double rank 0 items %d data follows'%(x_bins,y_bins,z_bins,x_edges[0]+delta/2.,y_edges[0]+delta/2.,z_edges[0]+delta/2.,delta,delta,delta,x_bins,y_bins,z_bins,nBins))    #,footer='object "density (all) [A^-3]" class field'
        with open(dx_file_name,'a') as W:
            if nBins%3 == 2:
                W.write('%.18e %.18e\nobject "density (all) [A^-3]" class field'%(append_able_data[0],append_able_data[1]))
            elif nBins%3 == 1:
                W.write('%.18e\nobject "density (all) [A^-3]" class field'%(append_able_data[0]))
        # deal with avg cv histogram
        if cv_data != []:
            print_able_data = flattened_cv[:-(nBins%3)]
            reshaped = np.reshape(print_able_data,(int(len(flattened_cv)/3),3))
            append_able_data = flattened_cv[-(nBins%3):]
            ## origin not shifted by half delta...
            #np.savetxt(cv_file_name,reshaped,comments='',header='object 1 class gridpositions counts %d %d %d\norigin %.18e %.18e %.18e\ndelta %.18e 0 0\ndelta 0 %.18e 0\ndelta 0 0 %.18e\nobject 2 class gridconnections counts %d %d %d\nobject 3 class array type double rank 0 items %d data follows'%(x_bins,y_bins,z_bins,x_edges[0],y_edges[0],z_edges[0],delta,delta,delta,x_bins,y_bins,z_bins,nBins))    #,footer='object "density (all) [A^-3]" class field'
            # origin shifted by half delta...
            np.savetxt(cv_file_name,reshaped,comments='',header='object 1 class gridpositions counts %d %d %d\norigin %.18e %.18e %.18e\ndelta %.18e 0 0\ndelta 0 %.18e 0\ndelta 0 0 %.18e\nobject 2 class gridconnections counts %d %d %d\nobject 3 class array type double rank 0 items %d data follows'%(x_bins,y_bins,z_bins,x_edges[0]+delta/2.,y_edges[0]+delta/2.,z_edges[0]+delta/2.,delta,delta,delta,x_bins,y_bins,z_bins,nBins))    #,footer='object "density (all) [A^-3]" class field'
            with open(cv_file_name,'a') as W:
                if nBins%3 == 2:
                    W.write('%.18e %.18e\nobject "density (all) [A^-3]" class field'%(append_able_data[0],append_able_data[1]))
                elif nBins%3 == 1:
                    W.write('%.18e\nobject "density (all) [A^-3]" class field'%(append_able_data[0]))
    else:
        # deal with space
        reshaped = np.reshape(weighted_flattened_space, (int(len(weighted_flattened_space)/3),3))
        ## origin not shifted by half delta...
        #np.savetxt(dx_file_name,reshaped,comments='',header='object 1 class gridpositions counts %d %d %d\norigin %.18e %.18e %.18e\ndelta %.18e 0 0\ndelta 0 %.18e 0\ndelta 0 0 %.18e\nobject 2 class gridconnections counts %d %d %d\nobject 3 class array type double rank 0 items %d data follows'%(x_bins,y_bins,z_bins,x_edges[0],y_edges[0],z_edges[0],delta,delta,delta,x_bins,y_bins,z_bins,nBins),footer='object "density (all) [A^-3]" class field')
        # origin shifted by half delta...
        np.savetxt(dx_file_name,reshaped,comments='',header='object 1 class gridpositions counts %d %d %d\norigin %.18e %.18e %.18e\ndelta %.18e 0 0\ndelta 0 %.18e 0\ndelta 0 0 %.18e\nobject 2 class gridconnections counts %d %d %d\nobject 3 class array type double rank 0 items %d data follows'%(x_bins,y_bins,z_bins,x_edges[0]+delta/2.,y_edges[0]+delta/2.,z_edges[0]+delta/2.,delta,delta,delta,x_bins,y_bins,z_bins,nBins),footer='object "density (all) [A^-3]" class field')
        # deal with cv
        if cv_data != []:
            reshaped = np.reshape(flattened_cv, (int(len(flattened_cv)/3),3))
            ## origin not shifted by half delta...
            #np.savetxt(cv_file_name,reshaped,comments='',header='object 1 class gridpositions counts %d %d %d\norigin %.18e %.18e %.18e\ndelta %.18e 0 0\ndelta 0 %.18e 0\ndelta 0 0 %.18e\nobject 2 class gridconnections counts %d %d %d\nobject 3 class array type double rank 0 items %d data follows'%(x_bins,y_bins,z_bins,x_edges[0],y_edges[0],z_edges[0],delta,delta,delta,x_bins,y_bins,z_bins,nBins),footer='object "density (all) [A^-3]" class field')
            # origin shifted by half delta...
            np.savetxt(cv_file_name,reshaped,comments='',header='object 1 class gridpositions counts %d %d %d\norigin %.18e %.18e %.18e\ndelta %.18e 0 0\ndelta 0 %.18e 0\ndelta 0 0 %.18e\nobject 2 class gridconnections counts %d %d %d\nobject 3 class array type double rank 0 items %d data follows'%(x_bins,y_bins,z_bins,x_edges[0]+delta/2.,y_edges[0]+delta/2.,z_edges[0]+delta/2.,delta,delta,delta,x_bins,y_bins,z_bins,nBins),footer='object "density (all) [A^-3]" class field')
    
    return half_max_counts

