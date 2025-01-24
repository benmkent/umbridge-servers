#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#"""

from subprocess import call
import os
import re

call(["rm","coords1"])
call(["rm","coords2"])
call(["rm","coords3"])
call(["rm","coords4"])
call(["rm","coords5"])

with open("../../coords/humpZone1.txt", 'r') as f1:
    lines = f1.readlines()
    for line in lines:
        line = "("+line.rstrip()+" 0)"+"\n"
        coordinates = re.sub('[;^M]',' ',line)
        with open("coords1", 'a') as outputFile:
            outputFile.write(coordinates)
        with open("coords10", 'a') as outputFile2:
            outputFile2.write(coordinates)
            
with open("../../coords/humpZone2.csv", 'r') as f2:
    lines = f2.readlines()
    for line in lines:
        line = "("+line.rstrip()+")"+"\n"
        coordinates = re.sub('[;^M]',' ',line)
        with open("coords2", 'a') as outputFile:
            outputFile.write(coordinates)
            
with open("../../coords/humpZone3.csv", 'r') as f3:
    lines = f3.readlines()
    for line in lines:
        line = "("+line.rstrip()+")"+"\n"
        coordinates = re.sub('[;^M]',' ',line)
        with open("coords3", 'a') as outputFile:
            outputFile.write(coordinates)
            
with open("../../coords/humpZone4.csv", 'r') as f4:
    lines = f4.readlines()
    for line in lines:
        line = "("+line.rstrip()+")"+"\n"
        coordinates = re.sub('[;^M]',' ',line)
        with open("coords4", 'a') as outputFile:
            outputFile.write(coordinates)

with open("../../coords/humpZone5.csv", 'r') as f5:
    lines = f5.readlines()
    for line in lines:
        line = "("+line.rstrip()+")"+"\n"
        coordinates = re.sub('[;^M]',' ',line)
        with open("coords5", 'a') as outputFile:
            outputFile.write(coordinates)
