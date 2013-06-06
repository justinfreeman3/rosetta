#Script Name: rosettaTranslate.py
#Description: This script is the second half of the Rosetta pipeline. It shifts bed file data from one genomic coordinate system to another.
#Author: Justin Freeman
#Organization: Dowell Laboratory, University of Colorado at Boulder
#Contact: justin.freeman@colorado.edu or robin.dowell@colorado.edu
#
#Version: 1.0, June 05 2013

#Imports
from __future__ import division
import os
import sys
import glob
import string
import argparse

parser = argparse.ArgumentParser(description='This script is the second half of the Rosetta pipeline.')
parser.add_argument('-g', '--masterGenomeFile', help='File containing cross-coordinate system genome data', required=True)
parser.add_argument('-o', '--outputFile', help='Desired filename prefix for output data', required=False)
parser.add_argument('-s', '--strand', help='Source strand for data, either plus or minus', required=True)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-q', '--queryFile', help='File containing data from the query genome (may not be used with -r)')
group.add_argument('-r', '--refFile', help='File containing data from the reference genome (may not be used with -q')

args=vars(parser.parse_args())

def main(args):
	print '\n\n' + '--------------------------------------------' + '\n'
	print 'Rosetta Translate is now running...' + '\n'	
	print '--------------------------------------------' + '\n'

	#Establish the output file prefix, if none given by the user
	if args['outputFile'] is None:
		args['outputFile'] = 'rosettaTranslate'
	else:
		pass 

	#Process the query genome data, if there is any
	if args['queryFile'] is not None:
		print 'Processing query genome data from: ' + args['queryFile'] + '\n'
		intermediateFile = processQueryData(args)
		print 'Query genome data processing complete. Stored in temporary file: ' + intermediateFile + '\n'
	else:
		pass
	
	#Process the reference genome data, if there is any
	if args['refFile'] is not None:
		print 'Processing reference genome data from: ' + args['refFile'] + '\n'
		intermediateFile = processReferenceData(args)
		print 'Reference genome data processing complete. Stored in temporary file: ' + intermediateFile + '\n'
	else:
		pass

	#Create a wig file from the intermediate files produced by either data processing function (query/reference) above
	print 'Creating a wig file from the processed data.' + '\n'	
	wigFile = makeWigFile(args, intermediateFile)
	
	#Convert the wig file into a bed file
	print 'Converting the wig file into a bed file.' + '\n'	
	bedFile = convertWigToBed(args, wigFile)

	print '\n' + '--------------------------------------------' + '\n'
	print 'Rosetta Translate is now complete!' + '\n'	
	print '--------------------------------------------' + '\n'


def processQueryData(args):
	queryDictionary = {}
	with open(args['queryFile']) as queryFile:
		for line in queryFile:
			(chromosome, position, value) = line.split()
			uniqueKey = chromosome + position
			queryDictionary[uniqueKey] = value
	intermediateFileName = args['outputFile'] + '.processQuery'
	intermediateFile = open(intermediateFileName, 'w')
	with open(args['masterGenomeFile']) as genomeFile:
		for line in genomeFile:
			genomeLine = line.split()
			qryChromosome = genomeLine[2]
			qryPosition = genomeLine[3]
			uniqueGenomeKey = qryChromosome + qryPosition
			if uniqueGenomeKey in queryDictionary:
				intermediateFile.write(line.rstrip() + '\t' + queryDictionary[uniqueGenomeKey] + '\n')
			else:
				intermediateFile.write(line.rstrip() + '\t' + '0' + '\n')
	intermediateFile.close()
	return intermediateFileName

def processReferenceData(args):
	print 'This is currently a dummy function. No work done here.' + '\n'
	intermediateFileName = 'test'
	return intermediateFileName

def makeWigFile(args, intermediateFile):
	currentChromosome = ''
	wigFileName = args['outputFile'] + '.wig'
	wigFile = open(wigFileName, 'w')
	with open(intermediateFile) as masterFile:
		for line in masterFile:
			dataLine = line.split()
			if dataLine[0] <> currentChromosome:
				currentChromosome = dataLine[0]
				wigFile.write('fixedStep chrom=' + currentChromosome + ' start=1 step=1\n')
			else:
				positionValue = dataLine[-1]
				wigFile.write(str(positionValue) + '\n')
	wigFile.close()
	return wigFileName

def convertWigToBed(args, wigFileName):
	bedFileName = args['outputFile'] + '.bed'
	bedFile = open(bedFileName, 'w')
	with open(wigFileName) as wigFile:
		regionLocation = 0
		inRegion = 0
		regionStart = 0
		regionEnd = 0
		strand = args['strand']
		if strand == 'plus':
			strand = '+'
			color = '238,0,0'
		else:
			strand = 'minus'
			color = '0,0,175'
		for line in wigFile:
			lineValue = line.split()
			if len(lineValue[0]) == 1:	#Determine if current line is a new chromosome header or not
				if int(lineValue[0]) == 0:	#Not in a feature region
					if inRegion == 1:	#Need to end the current feature region, print the output, and reset the feature region flag
						regionEnd = regionLocation
						printRegion(currentChromosome, regionStart, regionEnd, bedFile, strand, color)
						inRegion = 0
					else:	#Don't need to do anything, since we're already outside of a feature region
						pass
				elif int(lineValue[0]) == 1:	#In a feature region
					if inRegion == 0:	#Just now entering the feature region, need to set the feature region flag and the start location
						regionStart = regionLocation
						inRegion = 1
					else:	#Don't need to do anything, since we're just continuing with the same feature region
						pass
				regionLocation += 1		#Increment the location before moving to the next line of data in the wig file
			else:	#For a new chromosome
				if inRegion == 1:	#Need to wrap up and print the current feature region, since it clearly can't carry over to the next chromosome
					regionEnd = regionLocation
					printRegion(currentChromosome, regionStart, regionEnd, bedFile, strand, color)
					inRegion = 0
				else:	#Establish new chromosome value
					lineValue = line.split()
					currentChromosome = lineValue[1]
					currentChromosomeNumber = currentChromosome[6:]
					currentChromosome = currentChromosomeNumber.rstrip()
				regionLocation = 0
	bedFile.close()
	return bedFileName

def printRegion(chromosome, regionStart, regionEnd, bedFile, strand, color):
	bedFile.write(chromosome + '\t' + str(regionStart) + '\t' + str(regionEnd) + '\t' + '1\t0\t' + strand + '\t' + str(regionStart) + '\t' + str(regionEnd) + '\t' + color + '\n')


if __name__ == '__main__':
    main(args) 
