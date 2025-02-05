from functions import value_dict_maker, struct_dict_maker, lib_constructor
import sys
import os 

information = dict({
    "bases": ["def2-tzvp", "def2-svpd", "def2-svp"],
    "functionals": ["blyp", "pbe", "tpss"],
    "gcp": ["def2-tzvp", "def2-svpd", "def2-svp"]
})


# loop to acquire basis set, functional, and error function
while True:
    print("Enter appropriate values as instructed in documentation")
    basis = input("What basis set do you want: ").lower()
    functional = input("What functional do you want: ").lower()
    
    if basis in information["bases"] and functional in information["functionals"]:
        print(f"Basis set is {basis}. Functional is {functional}.\n")
        break
    else:
        print("There was a problem in entering the information. Check what you have input for correctness.\n")

# loop to acquire DFT only file, reference energy file and directory for xyz files
while True:
    print("Give all the following information as instructed")
    struct_path = input("Write out the full path for the directory that holds the structures' xyz files: ").lower()
    
    if os.path.isdir(struct_path):
        print("Directory confirmed to exist\n")
        struct_dict = struct_dict_maker(struct_path)

        if len(struct_dict) != 0:
            print(" All information matches is correct.\n")
            break

        print("Directory has no structures in it. Either rewrite path or select exit to correct directory")
        response = input("Do you want to rewrite the paths(y/n)? ").lower()
        if response == 'n' or response == 'no':
            sys.exit(1)
    else:
        print("The path provided does not lead to an existing directory - check paths for correctness\n")

gcp = False


if basis in information["gcp"]:
    correction = input("Do you want Geometric Counterpoise Corrections added(y/n): ").lower()
    


    if correction == "y" or "yes":
        gcp = True
else:
    print("No gCP available")

name = input("What name do you want to give the library: ")

lib_constructor(basis, functional, struct_dict, name, gcp)