from utils.MR import *
import glob

path2files='SMR.SCOPE13.CH01@Acquisition/' # path to MR data, they should be unziped (.sdds)

files_list = sorted(glob.glob(path2files+'*')) # create a list with all the MR data

for my_file in files_list: # create folder ./mr
    sdds_to_file(my_file)
    
make_pickle() # make pickle from all the sdds data (functions are in MR.py)