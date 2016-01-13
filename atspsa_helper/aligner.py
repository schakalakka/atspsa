import os

from config import *

dirs = os.listdir(DIR)

for filename in dirs:
    curr_DIR = DIR + filename
    os.chdir(curr_DIR)
    os.system("{}/daligner -l30 {}.db {}.db ".format(DALIGNER, filename, filename))
    os.system("{}/LAsort *.las".format(DALIGNER))
    os.system("{}/LAmerge {}.las *.S.las".format(DALIGNER, filename, filename))
    os.system("{}/LAdump -cdo {}.db {}.las > {}.ovl ".format(DALIGNER, filename, filename, filename))
