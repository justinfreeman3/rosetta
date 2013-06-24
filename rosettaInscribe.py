#Script Name: rosettaInscribe.py
#Description: This script is the first half of the Rosetta pipeline. It "inscribes," or produces the Rosetta file used to rapidly translate
#sequencing datasets between different genomic coordinate systems.
#Author: Justin Freeman
#Organization: Dowell Laboratory, University of Colorado at Boulder
#Contact: justin.freeman@colorado.edu or robin.dowell@colorado.edu
#
#Version: 1.0, June 11 2013

#Imports
from __future__ import division
import os
import sys
import glob
import string
import argparse
import operator

#-------------------------------------------------------------------------
#This section contains all the command line argument information, passed to the program by the user 
#and handled by the argparse module
#-------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='This script is the second half of the Rosetta pipeline.')


args=vars(parser.parse_args())
#-------------------------------------------------------------------------





#-------------------------------------------------------------------------

if __name__ == '__main__':
    main(args) 
