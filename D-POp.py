from functions import value_dict_maker, struct_dict_maker,energy, error_function, simplex, nelder_mead
from sklearn.metrics import mean_absolute_error as mae
from sklearn.metrics import mean_absolute_percentage_error as mape
from sklearn.metrics import root_mean_squared_error as rmse
import sys
import os 
import subprocess
from datetime import datetime
import ast

information = dict({
    "bases": ["def2-tzvp", "def2-svpd", "def2-svp"],
    "functionals": ["blyp", "pbe", "tpss"],
    "error_funcs": ["mae", "mape", "rmse", "set_error"]
})

# important initializing variables
param_num = 3

# condition to check for a "checkpoint" file - file designated to store state of simplex as it runs
if os.path.isfile("checkpoint"):
    restart = input("There is a checkpoint file which means a calculation did not finish. Would you like to continue the previous calculation (y/n): ").lower()
    
    if restart == "y" or restart == "yes":
        print("Ensure that you retype all the details (basis set, functional, error function, reference values and DFT energy files exactly as before")
        with open("checkpoint", 'r') as c:
            splx = ast.literal_eval([x for x in c.read().split("\n") if x][-1])
    else:
        splx = simplex(param_num)

    subprocess.run("rm checkpoint", shell=True)
else:   
    splx = simplex(param_num)
    

# loop to acquire basis set, functional, and error function
while True:
    print("Enter appropriate values as instructed in documentation")
    basis = input("What basis set do you want: ").lower()
    functional = input("What functional do you want: ").lower()
    error_func = input("What error function do you want: ").lower()
    if basis in information["bases"] and functional in information["functionals"] and error_func in information["error_funcs"]:
        print(f"Basis set is {basis}. Functional is {functional}. Error function is {error_func}\n")
        break
    else:
        print("There was a problem in entering the information. Check what you have input for correctness.\n")

# loop to acquire all reference files
while True:
    print("Give all the following information as instructed")
    lib_path = input("Write out the  full path to the file for the energies from the DFT only calculation: ").lower()
    ref_path = input("Write out the full path to the file for the reference energies: ").lower()
    struct_path = input("Write out the full path for the directory that holds the structures' xyz files: ").lower()

    # condition to check if reference files are present
    if os.path.isfile(lib_path) and os.path.isfile(ref_path) and os.path.isdir(struct_path):
        print("Files and directory confirmed to exist\n")
        lib_energies = value_dict_maker(lib_path)
        struct_dict = struct_dict_maker(struct_path)
        ref_energies = value_dict_maker(ref_path)

        # condition to check if reference files have the same number of structures
        if len(lib_energies) == len(ref_energies) == len(struct_dict):
            print(" All information matches in number of structures. Optimization shall begin soon\n")
            break

        print("Number of structures in references do not match. Either rewrite paths or select exit to correct files/directory")
        response = input("Do you want to rewrite the paths(y/n)? ").lower()
        if response == 'n' or response == 'no':
            sys.exit(1)
    else:
        print("One or both of the files do not exist or the path provided does not lead to an existing directory - check paths for correctness\n")

# for loop to remove any present orca input and junk files
for file in os.listdir(os.getcwd()):
    if "input" in file:
        subprocess.run("rm input*", shell=True)
        break


gcp = lib_path.split("_")[-1]

# initializing name of output file
file_name = f"output_{functional}_{basis}_{error_func}_{gcp}_{datetime.now().strftime("%m-%d-%Y_%H-%M-%S")}" 
file = open(file_name, 'x')
file.close()

with open(file_name, 'w') as f:
            f.write(f"{gcp}\n")

# running of Nelder-Mead Method
nelder_mead(splx, func=error_function, error_func=error_func, inter_energies= ref_energies, base_calculations=lib_energies, struct_dict=struct_dict, functional=functional, file=file_name)