# Python script to convert the output of ANGSD's -doGeno 5 function to dadi input. This script uses Python 2.7 and depends on having pandas and numpy installed. I usually edit this file in a text editor, save it, then run it in the command line by navigating to the directory where this script is found and typing: python convert_angsd_dogeno5_to_dadi.py

# Anything with ***** in the comment line requires user-specific values on that line

import pandas as pd
import numpy as np

# Point the the_file to your gunzipped -doGeno 5 file. You will have to add new columns for each of your individuals

geno = pd.read_csv("/Users/AirAlex/Documents/Berkeley/IBDIBE/reDone_analyses/genetics/sceloporus/redone_justPA_individuals/dadi/dogeno5_mccw_mingeno5.geno",delimiter="\t",index_col=False,names=["locus","position","ref","derived","ind1","ind2","ind3","ind4","ind5","ind6","ind7"]) # *****

# Create a blank data frame
new_data = []
# Make the column headers. Change the population names to ones relevant to your study (these must be the same ones that you use in subsequent dadi steps)
new_data.append(['Reference','Derived','Allele1','CW','MC',
                'Allele2','BD','MC','Gene','Position']) # *****


for row in geno.iterrows():
    data = row[1].tolist() #Turn each row into a list
    pop1 = " ".join(data[4:8]) # The first individual starts in column 4. Remember, in python, the last number in the range NOT included. So if your last individual for pop1 is in column 8, the range should be 4:9 (i.e. five individuals in pop1). Change these values to include all individuals in population 1 *****
    pop2 = " ".join(data[8:11]) # Change these values to include all individuals in population 2. The first number should be the same as the last number in pop1, and the second number should be one greater than the number of columns present in the doGeno 5 file. *****
    pop1_ref = pop1.count(data[2]) # Count the number of reference alleles (data[2] for this row) in pop 1
    pop2_ref = pop2.count(data[2])
    pop1_derived = pop1.count(data[3]) # Count the number of derived alleles (data[3]) in pop1
    pop2_derived = pop2.count(data[3])
    
    new_data.append(["-"+data[2]+"-", "-"+data[3]+"-",data[2],pop1_ref,pop2_ref,
                     data[3],pop1_derived,pop2_derived,data[0],data[1]]) # Put it all together


## For your convenience, this will output the genotypes of the last locus for both of the populations. Compare this to the doGeno 5 output to make sure you assigned your individuals correctly. It also prints one example line of the dadi input that you just created.
print pop1
print pop2
print new_data[1]

import csv

the_file = open("/Users/AirAlex/Documents/Berkeley/IBDIBE/reDone_analyses/genetics/sceloporus/redone_justPA_individuals/dadi/cwmc_dadi_input.csv", "w") # Change the output file to be whatever/wherever you would like *****
writer = csv.writer(the_file, delimiter = "\t") # Create the guts of the file
writer.writerows(new_data) # Put the data that you just made into the file
the_file.close() # Close the file
