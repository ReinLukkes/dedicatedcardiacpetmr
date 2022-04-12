# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 08:46:53 2020

@author: ischagt
"""
#import os
#from pathlib import Path
#file_name = open(r"C:\castor_v3.0.1\build\Release\umcupetmr\umcupetmr_translation_1_df.Cdf")

#number_of_lines = len(open(r"C:\castor_v3.0.1\build\Release\umcupetmr\umcupetmr_translation_1_df.Cdf").readlines(  ))
#print(number_of_lines)

#with open(rb"C:\castor_v3.0.1\build\Release\umcupetmr\umcupetmr_translation_1_df.Cdf") as file: # b is important -> binary
#    fileContent = file.read()

#open file in read mode
#file = open(r"C:\castor_v3.0.1\build\Release\umcupetmr\umcupetmr_translation_1_df.Cdf")

#read the content of file
#data = file.read()

#get the length of the data
#number_of_characters = len(data)

#print('Number of characters in text file :', number_of_characters)
#print(data[0,100,1])
#file_stats = os.stat(file_name)
#print(file_content)
#print(file_stats)
#print(f'File Size in Bytes is {file_stats.st_size}')
#lines=0
##for c in range(number_of_characters):
#    d = data.find("\n")
#    if d!=-1:
#        lines=lines+1
#string.find(value, start, end)
#print(data)
    
with open(rb"C:\castor_v3.0.1\build\Release\umcupetmr\umcupetmr_translation_1_df.Cdf") as f:
    byte = f.read(1)
    while byte != b"":
        print(byte)
        byte = f.read(1)