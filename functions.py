import subprocess
from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import mean_absolute_percentage_error as mape
from sklearn.metrics import root_mean_squared_error as rmse
import os
import pickle
from statistics import stdev
from random import random, uniform
import math


################################ INTERACTION ENERGY CALCULATION SPECFIC FUNCTIONS & VARIABLES ######################################################

#////// Variables ////////////////////////////////

# variables that anchor the paramter variables to be picked - see get params method
param_vars = [ "D3A1", "D3A2", "D3S8", "!"]

# 
# conversion factor for 1 hatree to kcals/mol
au = 627.5095 

# file that contains beginning of input file for OCRA
# this file is needed to construct the full input file 
# it contains information for the functional, dispersion parameters and type of dispersion correction
top = "/home/riley/Dispersion-Parameter-Optimization-D-POp-/tt"

# file that contains beginning of input file for OCRA
# this file is needed to construct the full input file 
# it contains information for the functional, basis set, and type of gcp for calculation of library for dispersion free calculations
lib_top = "/home/riley/Dispersion-Parameter-Optimization-D-POp-/top"

# file that contains end of input file for OCRA
# this file is needed to construct the full input file 
bottom = "/home/riley/Dispersion-Parameter-Optimization-D-POp-/bb"

# dictionary to hold possible gcp values
gcp_dict = dict({
    "def2-tzvp": "gcp(dft/tz)",
    "def2-svpd": "gcp(file)",
    "def2-svp" : "gcp(dft/svp)"
})
#////// Variables ////////////////////////////////

#////// Helper Functions ////////////////////////

# function to return list of parameters from the top of the orca input file
def get_params(file):
    file_list = file.split()
    return [file_list[file_list.index(y) + 1] for y in param_vars]


# function to create dictionary to hold structures
#struct_dicts is a dictionary {key: value} of {structure name: [dimer path, monomer 1 path, monomer 2 path]}
def struct_dict_maker(path):

    # initialization of 
    struct_dict = dict()
    file_list = [x for x in sorted(os.listdir(path)) if ".xyz" in x]
    
    for file in file_list:
        key = (file.replace(".xyz", "")
               .replace("dd", "")
               .replace("m1","")
               .replace("m2", "")
               .replace("a_","")
               .replace("c_", ""))
        if key in struct_dict.keys():
            struct_dict[key].append(f"{path}/{file}")
        else:
            struct_dict.update({key:[f"{path}/{file}"]})
            
    return struct_dict

# function to create dictionary to reference energy values
#dictionary {key: value} of {structure name: energy value}
def value_dict_maker(path):
    with open(path, "r") as f:
        value_list = [ x.split() for x in f.read().split("\n") if x]
    return {x[0]: float(x[1]) for x in value_list}

# function to change top of orca input file
def top_change(top_file, change_data):
    with open(top_file, "w") as t:
        t.write(change_data)

# function to acquire basis set, functional and gcp component of the top of the input file for library calculations -lib_top
def get_lib_params(file):
    lib_list = file.split()
    return [lib_list[lib_list.index(param_vars[-1]) + i] for i in [1, 2, 3]]

#////// Helper Functions ////////////////////////

# function to construct library of calculations that leave out dispersion correction - base_calculations in energy function
def lib_constructor(basis_set, functional, struct_dict, name, gcp):

    # opens up top of input file 
    with open(lib_top, "r") as t:
        file = t.read()

    # acquires paramters from top of input file and replaces them with inputs
    lib_params = get_lib_params(file)

    # checks if user wants gcp
    if gcp:
        change_file = (file.replace(lib_params[0], functional)
                    .replace(lib_params[1], basis_set)
                    .replace(lib_params[2], gcp_dict[basis_set]))
    else:
        change_file = (file.replace(lib_params[0], functional)
                    .replace(lib_params[1], basis_set)
                    .replace(lib_params[2], ""))        

    # changes file in directory t contain given parameters 
    top_change(top, change_file)

    # outer loop goes over all the keys - structure names 
    for struct in struct_dict.keys():

        # list to store DFT energies for dimer and monomers
        energies = []
        
        # inner loop goes over all items in list that contains paths for given structure name
        for path in struct_dict[struct]:

            # checks if name of structure has markers for anion or cation to change a specific value in input file
            if "a_" in path.split("/")[-1]:
                change_data = change_data.replace("xyz 0", "xyz -1")
            elif "c_" in path.split("/")[-1]:
                change_data = change_data.replace("xyz 0", "xyz 1")

            # constructs entire input file
            subprocess.run(f"cat {top} {path} {bottom} > input", shell=True)

            # runs energy calculation on ORCA, picks out dispersion energy and appends it to list
            energies.append(float(subprocess.check_output(f"/home/riley/bin/orca_6_0_0_shared_openmpi416/orca input| grep 'FINAL'| awk '{{print $5}}' ORS=' '",
                                                              shell=True,
                                                              executable="/bin/bash")
                                                              .decode()))

            # clears out all junk files for ORCA calculation
            subprocess.run("rm input*", shell=True)

        # opens file to append new entry ito library of energies
        with open(f"{name}_base_calculations_{functional}_{basis_set}_gcp={gcp}".lower(), 'a') as lib:
            lib.write(f"{struct} {au*(energies[0] - (energies[1]+energies[2]))}\n")
            
    top_change(lib_top, file)


# function to calculate the interaction energies of the dimers in the dataset given
def energy(change_params, kwargs):

    # necessary variables picked from the keyword argument (kwargs) dictionary
    struct_dict = kwargs["struct_dict"]
    inter_energies = kwargs["inter_energies"]
    base_calculations = kwargs["base_calculations"]
    functional = kwargs["functional"]
    
    # opens up top file 
    with open(f"{top}", "r") as f:
        data = f.read()

    # change_params contains the dispersion parameters we want to test
    # this portion converts them to strings and adds on the functional so they can be replaced at once
    change = [str(y) for y in change_params]
    change.append(functional)

    # this portion acquires the original paramters [including the functional] 
    # the original parameters are changed for our input parameters
    params = get_params(data)
    change_data = (data.replace(params[0], change[0])
                    .replace(params[1], change[1])
                    .replace(params[2], change[2])
                    .replace(params[3], change[3]))

    # these are lists to store the interaction enegies calculated 
    # and the reference energies that were found with more accurate methods
    energy_list = []
    strd_values = []

    
    # outer loop cycles through all of the keys - structure names 
    for struct in struct_dict.keys():

        # list to store dispersion energies for dimer and monomers
        energies = []

        # appending to reference value list
        strd_values.append(inter_energies[struct])

        # inner loop goes over all items in list that contains paths for given structure name
        for path in struct_dict[struct]:

            # checks if name of structure has markers for anion or cation to change a specific value in input file
            if "a_" in path.split("/")[-1]:
                change_data = change_data.replace("xyz 0", "xyz -1")
            elif "c_" in path.split("/")[-1]:
                change_data = change_data.replace("xyz 0", "xyz 1")

            # changes top of input file to contain appropriate values 
            top_change(top, change_data)

            # constructs entire input file
            subprocess.run(f"cat {top} {path} {bottom} > input", shell=True)

            # runs energy calculation on ORCA, picks out dispersion energy and appends it to list
            energies.append(float(subprocess.check_output(f"/home/riley/bin/orca_6_0_0_shared_openmpi416/orca input| grep 'Dispersion'| awk '{{print $3}}' ORS=' '",
                                                              shell=True,
                                                              executable="/bin/bash")
                                                              .decode()))

            # clears out all junk files for ORCA calculation
            subprocess.run("rm input*", shell=True)

            # changes top of input file back to its original state 
            # to prevent information from bleeding between calculations
            top_change(top, data)

        # calculates interaction energy and appends it to energy list
        energy_list.append(base_calculations[struct] + au*(energies[0] - (energies[1]+energies[2])))

    return strd_values, energy_list

################################ INTERACTION ENERGY CALCULATION SPECFIC FUNCTIONS & VARIABLES ######################################################


################################ NELDER MEAD METHOD SPECFIC FUNCTIONS & VARIABLES ##################################################################

#////// Helper Functions ////////////////////////

# function to sort two lists by the values in the first 
def sorter(values, splx):
    y, x = zip(*sorted(zip(values,splx)))
    return list(x), list(y)

# function to create random simplex with points in the range of 0-5 in stated dimensions
def simplex(dim):
    splx = []
    for i in range(dim+1):
        splx.append([(x/x)*uniform(0, 5) for x in range(1,dim+1)])
    return splx

# function to calculate the centroid/mean of vertices of simplex except that with the highest value
def mean(splx):
    return [sum(y)/len(y) for y in list(zip(*splx))]

#////// Helper Functions ////////////////////////

#////// Error Functions ////////////////////////

# general error function
def error_function(change_params, kwargs):
    strd_values, energy_list = energy(change_params,kwargs)
    match kwargs["error_func"]:
        case "mae":
            return mae(strd_values, energy_list)

        case "mape":
            return mape(strd_values, energy_list)

        case "rmse":
            return rmse(strd_values, energy_list)

        case "set_error":
            error = 1
            deviations = [abs(x - y) for (x,y) in zip(strd_values, energy_list)]
            for dv in deviations:
                error = error * dv
            return error

#////// Error Functions ////////////////////////

#Nelder-Mead Optimization - referenced from Algorithms for Optimization by Mykel J. Kochenderfer and Tim A. Wheeler
def nelder_mead( splx, α=1.0, β=2.0,γ=0.5, ϵ=0.0001, **kwargs):

    # function to minimize in the algorithm 
    func = kwargs["func"]
    
    # number to keep track of what iteration of the following loop is currently running
    count = 0
    
    # Δ is the convergence condition variable - in this case it is the standard deviation of the values 
    # on the error function
    Δ = float('inf')
    print("Calculating error values of simplex points")

    # calculates values for each point on the simplex from the error function
    y_vals = [func(x, kwargs) for x in splx]

    # loop runs until standard deviation is lower than the given tolerance ϵ
    while Δ > ϵ:

        # sorts the simplex vertex points based on values in error function
        splx, y_vals = sorter(y_vals,splx)

        # updates count of iterations 
        count = count + 1

        # updates checkpoint file with latest simplex
        with open("checkpoint", 'a') as c:
            c.write(f"Iteration {count}\n{str(splx)}\n\n")
            
        # extracts lowest point and its value
        xl, yl = splx[0], y_vals[0]

        # extracts second highest point and its value
        xs, ys = splx[-2], y_vals[-2]

        # extracts highest point and its value
        xh, yh = splx[-1], y_vals[-1]
        print("Calculating mean point")

        # calculates centroid/mean of all points except the highest
        xm = mean(splx[0:len(splx)-1])

        # calculates reflection of highest point across the centroid/mean % its value
        xr = [a + b for (a,b) in zip(xm,[α*(c-d) for (c,d) in zip(xm,xh)])] #xm + (α*(xm-xh))
        yr = func(xr, kwargs)

        # condition to check if reflected point value is lower than lowest value
        if yr < yl:
            print("Testing extension point")

            # if the conditon is passed it means this direction is favorable
            # an extension point is found further in the direction of the reflected point
            xe = [a + b for (a,b) in zip(xm,[β*(c-d) for (c,d) in zip(xr,xm)])] #xm + (β*(xr-xm))
            ye = func(xe, kwargs)

            # if the value of the extension point is better than the reflected point, it replaces the highest point
            # if it is higher the reflected point, the reflected point replaces the highest point
            # this means the formerly second highest point becomes the highest point in the next iteration
            splx[-1], y_vals[-1] = (xe,ye) if ye < yr else (xr,yr)

        # outer condition to check if reflected point value is greater than or equal to second highest point
        elif yr >= ys:
            print("Testing contraction point")

            # if condition is met, this condtion checks if it is lower than the highest point
            if yr < yh:
                # replaces highest point if condition is true
                xh, yh, splx[-1], y_vals[-1] = xr, yr, xr, yr

            # if outer condition is met, it implies that going further out is not worthwhile
            # a contraction point and its value is calculated 
            xc = [a + b for (a,b) in zip(xm,[γ*(c-d) for (c,d) in zip(xh,xm)])] #xm + (γ*(xh-xm))
            yc = func(xc, kwargs)

            # condition to check if contraction point value is larger than highest value
            if yc > yh:
                print("Squeezing simplex to smaller size")

                # if condition is met it means contracting inward towards the simplex is not favorable
                # the only option left is to shrink the simplex towards the lowest point
                for i in range(1, len(y_vals)):
                    splx[i] = [(a + b)/2 for (a,b) in zip(splx[i],xl)] #(splx[i] + xl)/2
                    y_vals[i] = func(splx[i], kwargs)
            else:
               # if condition is not met the highest point is replaced by the contraction point
               splx[-1], y_vals[-1] = xc, yc
        else:
            # if all condtions are not met simply replace the highest point with the reflected point
            # this effectively flips over the simplex and allows it to roam around
            splx[-1], y_vals[-1] = xr, yr

        # upating of standard deviation of point values
        Δ = stdev(y_vals)
        
    print("Finished")

    # command to delete checkpoint file
    subprocess.run("rm checkpoint", shell=True)
    
    # when loop is broekn the lowest point on the simplex is given as the parameters and the associated function value is the error
    error = min(y_vals)
    parameters = splx[y_vals.index(error)]

    # variables needed to print to output file
    strd_values, energy_list = energy(parameters, kwargs)
    file_name = kwargs["file"]
    error_func = kwargs["error_func"]

    # putting final information in output file
    with open(file_name, 'a') as f:
        f.write(f"The MAE is {mae(strd_values, energy_list)}. The MAPE is {mape(strd_values, energy_list)}. The RMSE is {rmse(strd_values, energy_list)}\nThe final parameters are\n A1: {parameters[0]}\n S8: {parameters[1]}\n A2: {parameters[2]}\n The {error_func.upper()} with these parameters is {error}\nThe list of interaction energies with the final parameters is:\n {energy_list}\n\n")
        
################################ NELDER MEAD METHOD SPECFIC FUNCTIONS & VARIABLES ##################################################################