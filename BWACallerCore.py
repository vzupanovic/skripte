#!/usr/bin/python
import subprocess
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import analyzer as stats
import operator
import matplotlib.ticker as mticker
from pylab import *
from BWACaller import *

class CallerCore:
	def __init__(self, readType, readFile_1, readFile_2, contigFile, plotPath, showPlot):
		self.maxCovReg = [-1,-1,-1,-1]
		self.readType = readType
		self.readFile_1 = readFile_1
		self.readFile_2 = readFile_2
		self.contigFile = contigFile
		self.plotPath = plotPath
		self.showPlot = showPlot
		
	def doBWA(self, regularPath):
		self.bwa = BWACaller(self.readType, self.contigFile, self.readFile_1, self.readFile_2, regularPath)
		self.bwa.calculateIndex()
		self.bwa.align()
		if self.readType == 1:
			self.bwa.doSamsa()
		else:
			self.bwa.doSampe()
		
	def doSAM(self, regularPath):
		print "[AN:] Getting basic alignment stats..."
		self.samtools = SAMTools(regularPath)
		self.samtools.execute()
		
	def doBAM(self, regularPath, pathToBedTools):
		self.bedtools = BEDTools(regularPath, pathToBedTools) #change this !
		self.bedtools.compute()
		
	def doCoverageStats(self, regularPath):
		print "[AN:] Getting basic coverage stats..."
		coverage_stats = Coverage()
		coverage_stats.getContigData(self.contigFile)
		coverage_stats.getAlignmentData(regularPath+"test.bam.cov")
		coverage_stats.getCoveragePerContig()
		potCollaps = coverage_stats.getPotColapseRegions(4)
		self.maxCovReg = potCollaps[:]
		print "[AN:] Potential colapse regions because of high coverage:"
		for contig in potCollaps:
			print "\tContig ID:",contig[0],"\tcoverage:",contig[1]
		print "[AN:] Plotting coverage per contig..."
		coverage_stats.plotAllContigCov(self.plotPath, self.showPlot)
		print "[AN:] Done..."
		
	def getAllStats(self):
		return self.maxCovReg
		

		

if __name__ == "__main__":
	if len(sys.argv) < 9:
		print "[AN:] Usage:\n\tpython BWACallerCore.py read1.fastq read2.fastq[null]\n\tcontig_file.fasta plot_folder read_type[0/1] show_plot[0/1]"
		print "\tpath_to_data_folder[~/example/] path_to_BEDTOOLS[~/bedtools/bin/bedtools]" 
	else:
		if os.path.isfile(sys.argv[1])!=True:
			print "[AN:] File "+str(sys.argv[1])+" doesn't exist. Exiting..."
			exit(-1)
		if os.path.isfile(sys.argv[2])!=True:
			print "[AN:] File "+str(sys.argv[2])+" doesn't exist. Exiting..."
			exit(-1)
		if os.path.isfile(sys.argv[3])!=True:
			print "[AN:] File "+str(sys.argv[1])+" doesn't exist. Exiting..."
			exit(-1)

		read1 = sys.argv[1]
		read2 = sys.argv[2]
		contigFile = sys.argv[3]
		plotFolder = sys.argv[4]
		readType = int(sys.argv[5])	
		showPlot = int(sys.argv[6])
		regularPath = sys.argv[7]
		pathToBedTools = sys.argv[8]
		#regularPath = "/home/lexy/BWA_test/"
		core = CallerCore(readType, read1, read2, contigFile, plotFolder, showPlot)
		core.doBWA(regularPath)
		core.doSAM(regularPath)
		core.doBAM(regularPath, pathToBedTools)
		core.doCoverageStats(regularPath)
		print core.getAllStats()
