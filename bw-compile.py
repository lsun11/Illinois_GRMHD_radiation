#!/bin/python
############################################################
# read BW compilation stdout and reorganize linker command #
# EXAMPLE: #################################################
# make -j8 giggle 1> compilation.out 2> compilation.err    #
# python ~/py/bw-compile.py compilation.out                #
############################################################

import sys,os,commands

## get the original (wrong) linker command from stdout file

FILE = sys.argv[1]
fd = open(FILE,'r')
LINES = fd.readlines()
fd.close()
for LINE in LINES:
    if LINE[:2]=="-L":
        ORIGINAL_LINKER_COMMAND = LINE.strip('\n')

#print ORIGINAL_LINKER_COMMAND
NEW_LINKER_COMMAND = "CC -o"+ORIGINAL_LINKER_COMMAND.split(' -o ')[1]+" "+ORIGINAL_LINKER_COMMAND.split(' -o ')[0]
print "================================================\n",\
       NEW_LINKER_COMMAND,\
      "\n================================================\n"

os.system(NEW_LINKER_COMMAND)

print "================================================\n"

OUTPUT = commands.getoutput("ls -l exe/cactus_giggle")

print OUTPUT
