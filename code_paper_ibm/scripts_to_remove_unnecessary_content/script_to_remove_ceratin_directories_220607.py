# -*- coding: utf-8 -*-
"""
Script to remove some unnecesary files from results to save space

"""

#REMOVE DIRECTORIES
import shutil
import os

path = "/Users/pablo/Desktop/code_paper/data/210822/"



for root, dirs, files in os.walk(path):
    for d in dirs:
        if  d == "per_round":
            #print("ok")
            filepath = root + os.sep + d
            filepath=filepath.replace("\\" , "/")
            #this is needed to remove a directory with content inside
            shutil.rmtree(filepath)
            

#REMOVE FILES
import os


path = "/Users/pablo/Desktop/code_paper/data/210822/"

for root, dirs, files in os.walk(path):
    for file in files:
        if file == "df_mutants.csv":
            #print("ok")
            filepath = root + os.sep + file
            filepath=filepath.replace("\\" , "/")
            os.remove (filepath)
        elif "grid0_round_" in file: #remove all the grid0_except the last one
            if "_49" in file:
                pass
            else:
                filepath = root + os.sep + file
                filepath=filepath.replace("\\" , "/")
                os.remove (filepath)
       

            