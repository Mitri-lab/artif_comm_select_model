# -*- coding: utf-8 -*-
"""
Script to remove some unnecesary files from results to save space

"""

#REMOVE DIRECTORIES
import shutil
import os

path = "/Users/pablo/Desktop/code_paper/data/210823_invasion/"



for root, dirs, files in os.walk(path):
    for d in dirs:
        if d == "d_s" or d=="d_r" or d=="d2_s" or d == "d2_r" or d == "0_25":
            filepath = root + os.sep + d
            filepath=filepath.replace("\\" , "/")
            #this is needed to remove a directory with content inside
            shutil.rmtree(filepath)
        elif  d == "per_round":
            #print("ok")
            filepath = root + os.sep + d
            filepath=filepath.replace("\\" , "/")
            #this is needed to remove a directory with content inside
            shutil.rmtree(filepath)
            

#REMOVE FILES
import os


path = "/Users/pablo/Desktop/code_paper/data/210823_invasion/"

for root, dirs, files in os.walk(path):
    for file in files:
        if "df_grid" in file:
            continue
        elif "species" in file:
            continue
        elif "INFO" in file:
            continue
        elif "index_rep" in file:
            continue
        elif "README" in file:
            continue
        elif "model" in file:
            continue
        else:
            filepath = root + os.sep + file
            filepath=filepath.replace("\\" , "/")
            os.remove (filepath)
        
   
       
   
       

            