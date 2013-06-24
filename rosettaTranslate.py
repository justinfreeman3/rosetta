#Script Name: rosettaTranslate.py
#Description: This script is the second half of the Rosetta pipeline. It shifts bed file data from one genomic coordinate system to another.
#Author: Justin Freeman
#Organization: Dowell Laboratory, University of Colorado at Boulder
#Contact: justin.freeman@colorado.edu or robin.dowell@colorado.edu
#
#Version: 1.2, June 11 2013

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
parser.add_argument('-g', '--masterGenomeFile', help='File containing cross-coordinate system genome data', required=True)
parser.add_argument('-o', '--outputFile', help='Desired filename prefix for output data', required=False)
parser.add_argument('-f', '--outputFormat', help='Desired file output format, (bed or bedgraph)', required=True)
#parser.add_argument('-s', '--strand', help='Source strand for data, either plus or minus', required=True)

#group = parser.add_mutually_exclusive_group(required=True)
parser.add_argument('-p', '--plusQueryFile', help='File containing plus strand data from the query genome (may not be used with -r)')
parser.add_argument('-m', '--minusQueryFile', help='File containing minus strand data from the query genome (may not be used with -r')
#group.add_argument('-r', '--refFile', help='File containing data from the reference genome (may not be used with -q')

args=vars(parser.parse_args())
#-------------------------------------------------------------------------

def main(args):
	#-------------------------------------------------------------------------
	#This function models the workflow for the program. No processing is done here. Functions are 
	#called as needed to provide processing and file output functionality
	#-------------------------------------------------------------------------

	print '\n\n' + '--------------------------------------------' + '\n'
	print 'Rosetta Translate is now running...' + '\n'	
	print '--------------------------------------------' + '\n'

	#Establish the output file prefix, if none given by the user
	if args['outputFile'] is None:
		args['outputFile'] = 'rosettaTranslate'
	else:
		pass 

	#Process the query genome data, if there is any
	if args['plusQueryFile'] is not None:
		print 'Processing query genome data from: ' + args['plusQueryFile'] + ' and ' + args['minusQueryFile'] + '\n'
		intermediateFile = processQueryData(args)
		print 'Query genome data processing complete. Stored in temporary file: ' + intermediateFile + '\n'
	else:
		pass
	
	#Process the reference genome data, if there is any
	#if args['refFile'] is not None:
	#	print 'Processing reference genome data from: ' + args['refFile'] + '\n'
	#	intermediateFile = processReferenceData(args)
	#	print 'Reference genome data processing complete. Stored in temporary file: ' + intermediateFile + '\n'
	#else:
	#	pass

	#Create a wig file from the intermediate files produced by either data processing function (query/reference) above
	print 'Creating a wig file from the processed data.' + '\n'	
	wigPlusFile, wigMinusFile = makeWigFile(args, intermediateFile)
	
	#Convert the wig file into a bed file or bedgraph file, based on input
	outputFormat = args['outputFormat']
	if outputFormat == 'bed':
		print 'Converting the wig file into a bed file.' + '\n'	
		bedPlusFile, bedMinusFile = convertWigToBed(args, wigPlusFile, wigMinusFile)
	else:
		print 'Converting the wig file into a bedgraph file.' + '\n'
		bedgraphPlusFile, bedgraphMinusFile = convertWigToBedgraph(args, wigPlusFile, wigMinusFile)

	print '\n' + '--------------------------------------------' + '\n'
	print 'Rosetta Translate is now complete!' + '\n'	
	print '--------------------------------------------' + '\n'
#-------------------------------------------------------------------------

def processQueryData(args):
	#-------------------------------------------------------------------------
	#This function ports data aligned to a query genome to the reference genome's 
	#coordinate system. It does this by creating a dictionary, wherein each entry's key
	#equals the concatenated string of 'chromosome' and 'position'. The value associated 
	#with each key is the actual experimental data value for that position. After generating
	#this dictionary, the program iterates through the rosetta file, searching for matches in 
	#the dictionary. If a match is found, the experimental data value is written into the 
	#rosetta file. If not, a value of '0' is written for that position.
	#-------------------------------------------------------------------------
		
	queryPlusDictionary = {}	#create dictionary to hold plus strand data
	queryMinusDictionary = {}	#create dictionary to hold minus strand data
	with open(args['plusQueryFile']) as plusQueryFile:
		for line in plusQueryFile:	#read plus strand data into dictionary
			lineValue = line.split()
			chromosome = lineValue[0]
			position = lineValue[1]
			value = lineValue[2]
			uniqueKey = chromosome + '_' + position
			queryPlusDictionary[uniqueKey] = value
	with open(args['minusQueryFile']) as minusQueryFile:	
		for line in minusQueryFile:	#read minus strand data into dictionary
			lineValue = line.split()
			chromosome = lineValue[0]
			position = lineValue[1]
			value = lineValue[2]
			uniqueKey = chromosome + '_' + position
			queryMinusDictionary[uniqueKey] = value

	intermediateFileName = args['outputFile'] + '.processQuery'
	intermediateFile = open(intermediateFileName, 'w')
	with open(args['masterGenomeFile']) as genomeFile:
		for line in genomeFile:
			genomeLine = line.split()
			qryChromosome = genomeLine[2]
			qryPosition = genomeLine[3]
			lastChromosome = ''
			lastPosition = ''			
			uniqueGenomeKey = qryChromosome + '_' + qryPosition
			if uniqueGenomeKey in queryPlusDictionary:
				plusValue = queryPlusDictionary[uniqueGenomeKey]
			else:
				plusValue = '0'
			if uniqueGenomeKey in queryMinusDictionary:				
				minusValue = queryMinusDictionary[uniqueGenomeKey]
			else:
				minusValue = '0'
			if qryChromosome == lastChromosome:
				inversionCheck = int(qryPosition) - int(lastPosition)				
				if inversionCheck < 0:	#This position is in an inversion					
					intermediateFile.write(line.rstrip() + '\t' + str(minusValue) + '\t' + str(plusValue) + '\n')						
				else:	#This position is not in an inversion				
					intermediateFile.write(line.rstrip() + '\t' + str(plusValue) + '\t' + str(minusValue) + '\n')
					#foundRead = 'FOUND' + uniqueGenomeKey
					#queryDictionary[foundRead] = queryDictionary.pop(uniqueGenomeKey)			
			else:
				intermediateFile.write(line.rstrip() + '\t' + str(plusValue) + '\t' + str(minusValue) + '\n')
			lastChromosome = qryChromosome
			lastPosition = qryPosition
	intermediateFile.close()
	#handleDanglingReads(args, queryDictionary)
	return intermediateFileName
#-------------------------------------------------------------------------

def processReferenceData(args):
	print 'This is currently a dummy function. No work done here.' + '\n'
	intermediateFileName = 'test'
	return intermediateFileName
#-------------------------------------------------------------------------

def handleDanglingReads(args, queryDictionary):
	print 'Trying to find dangling reads.' +'\n'
	numberOfDanglingReads = 0
	danglingDictionary = {}	#Not used anywhere yet
	summaryDictionary = {} 
	danglingFileName = args['outputFile'] + '_unmapped'
	danglingFile = open(danglingFileName, 'w')
	with open(args['queryFile']) as queryFile:
		for key in sorted(queryDictionary.iterkeys()):					
			if key[0:5] <> 'FOUND':
				spacer = key.find('_')
				chromosome = key[0:spacer]
				position = key[(spacer+1):]				
				danglingFile.write(chromosome + '\t' + position + '\t' + str(queryDictionary[key]) + '\n')
				numberOfDanglingReads += 1
				danglingDictionary
				if chromosome in summaryDictionary:
					summaryDictionary[chromosome] += 1
				else:
					summaryDictionary[chromosome] = 1
			else:
				pass
		print str(numberOfDanglingReads) + ' dangling reads found.' + '\n'
		#print summaryDictionary
	danglingFile.close()				
	return

#-------------------------------------------------------------------------

def makeWigFile(args, intermediateFile):
	currentChromosome = ''
	wigPlusFileName = args['outputFile'] + '_plus.wig'
	wigMinusFileName = args['outputFile'] + '_minus.wig'
	wigPlusFile = open(wigPlusFileName, 'w')
	wigMinusFile = open(wigMinusFileName, 'w')
	with open(intermediateFile) as masterFile:
		for line in masterFile:
			dataLine = line.split()
			if dataLine[0] <> currentChromosome:
				currentChromosome = dataLine[0]
				wigPlusFile.write('fixedStep chrom=' + currentChromosome + ' start=1 step=1' + '\n')
				wigMinusFile.write('fixedStep chrom=' + currentChromosome + ' start=1 step=1' + '\n')
			else:
				plusValue = dataLine[-2]
				minusValue = dataLine[-1]
				wigPlusFile.write(str(plusValue) + '\n')
				wigMinusFile.write(str(minusValue) + '\n')
	wigPlusFile.close()
	wigMinusFile.close()	
	return wigPlusFileName, wigMinusFileName
#-------------------------------------------------------------------------

def convertWigToBed(args, wigPlusFileName, wigMinusFileName):
	bedPlusFileName = args['outputFile'] + '_plus.bed'
	bedFile = open(bedPlusFileName, 'w')	
	wigFileName = wigPlusFileName
	
	counter = 0
	while counter < 2:	
		with open(wigFileName) as wigFile:
			regionLocation = 1
			inRegion = 0
			regionStart = 0
			regionEnd = 0
			if counter == 0:
				strand = '+'
				color = '238,0,0'
			else:
				strand = '-'
				color = '0,0,175'			
			for line in wigFile:
				lineValue = line.split()			
				if len(lineValue) == 1:	#Determine if current line is a new chromosome header or not
					if float(lineValue[0]) == 0:	#Not in a feature region
						if inRegion == 1:	#Need to end the current feature region, print the output, and reset the feature region flag
							regionEnd = regionLocation
							printRegion(currentChromosome, regionStart, regionEnd, bedFile, strand, color, regionValue)
							inRegion = 0
						else:	#Don't need to do anything, since we're already outside of a feature region
							pass
					elif float(lineValue[0]) <> 0:	#In a feature region
						if inRegion == 0:	#Just now entering the feature region, need to set the feature region flag and the start location
							regionStart = regionLocation + 1
							regionValue = lineValue[0]
							inRegion = 1
						else:	
							if regionValue == lineValue[0]:		#Don't need to do anything, since we're just continuing with the same feature region
								pass
							else:
								regionEnd = regionLocation
								printRegion(currentChromosome, regionStart, regionEnd, bedFile, strand, color, regionValue)
								inRegion = 0
					regionLocation += 1		#Increment the location before moving to the next line of data in the wig file
				else:	#For a new chromosome
					if inRegion == 1:	#Need to wrap up and print the current feature region, since it clearly can't carry over to the next chromosome
						regionEnd = regionLocation
						printRegion(currentChromosome, regionStart, regionEnd, bedFile, strand, color, regionValue)
						inRegion = 0
					else:	#Establish new chromosome value					
						currentChromosome = lineValue[1]
						currentChromosomeNumber = currentChromosome[6:]
						currentChromosome = currentChromosomeNumber.rstrip()
					regionLocation = 1		
		bedFile.close()
		if counter < 1:
			wigFileName = wigMinusFileName	
			bedMinusFileName = args['outputFile'] + '_minus.bed'
			bedFile = open(bedMinusFileName, 'w')
		else:
			pass
		counter += 1
	bedFile.close()
	return bedPlusFileName, bedMinusFileName
#-------------------------------------------------------------------------

def convertWigToBedgraph(args, wigPlusFileName, wigMinusFileName):
	bedgraphPlusFileName = args['outputFile'] + '_plus.bedgraph'
	bedgraphFile = open(bedgraphPlusFileName, 'w')
	bedgraphFile.write('type=bedGraph' + '\n')	#Print bedgraph file header
	wigFileName = wigPlusFileName
	counter = 0
	while counter <2:	
		with open(wigFileName) as wigFile:
			if counter == 0:
				strand = 'plus'
			else:
				strand = 'minus'
			for line in wigFile:
				lineValue = line.split()
				if len(lineValue) > 1:	#Establish new chromosome value
					currentChromosome = lineValue[1]
					currentChromosomeNumber = currentChromosome[6:]
					currentChromosome = currentChromosomeNumber.rstrip()
					regionStart = 0
					currentValue = ''
					currentPosition = 0
					lastValue = ''
				else:
					if lastValue == '':
						lastValue = lineValue[0]
						currentValue = lineValue[0]
					else:
						pass				
					if lastValue == lineValue[0]:	#Continue existing region
						pass
					else:	#Define new region
						regionEnd = currentPosition
						printBedgraphRegion(args, currentChromosome, regionStart, regionEnd, lastValue, bedgraphFile, strand)
						regionStart = currentPosition
						currentValue = lineValue[0]
				currentPosition += 1
				lastValue = currentValue	
		
		bedgraphFile.close()
		if counter < 1:
			wigFileName = wigMinusFileName
			bedgraphMinusFileName = args['outputFile'] + '_minus.bedgraph'
			bedgraphFile = open(bedgraphMinusFileName, 'w')
		counter += 1
	bedgraphFile.close()
	return bedgraphPlusFileName, bedgraphMinusFileName

#type=bedGraph
#chr1    0       4302    0
#chr1    4302    4303    0.0539228095498952
#chr1    4303    4304    0.0519969949231132
#-------------------------------------------------------------------------

def printRegion(chromosome, regionStart, regionEnd, bedFile, strand, color, regionValue):
	bedFile.write(chromosome + '\t' + str(regionStart) + '\t' + str(regionEnd) + '\t' + str(regionValue) + '\t' + '0' + '\t' + strand + '\t' + str(regionStart) + '\t' + str(regionEnd) + '\t' + color + '\n')
#-------------------------------------------------------------------------

def printBedgraphRegion(args, currentChromosome, regionStart, regionEnd, lastValue, bedgraphFile, strand):
	if strand == 'minus':	#Set value negative for minus strand regions
		if lastValue <> '0':
			if lastValue[0] == '-':
				pass
			else:			
				lastValue = '-' + lastValue
		else:
			pass
	else:
		pass
	bedgraphFile.write(currentChromosome + '\t' + str(regionStart) + '\t' + str(regionEnd) + '\t' + str(lastValue) + '\n')
#-------------------------------------------------------------------------

if __name__ == '__main__':
    main(args) 
