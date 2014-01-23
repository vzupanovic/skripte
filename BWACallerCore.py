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
	def __init__(self, sysArgs):
		if os.path.isfile(sysArgs[1])!=True:
			print "[AN:] File "+str(sysArgs[1])+" doesn't exist. Exiting..."
			exit(-1)
		if os.path.isfile(sysArgs[2])!=True:
			print "[AN:] File "+str(sysArgs[2])+" doesn't exist. Exiting..."
			exit(-1)
		if os.path.isfile(sysArgs[3])!=True:
			print "[AN:] File "+str(sysArgs[1])+" doesn't exist. Exiting..."
			exit(-1)
		self.maxCovReg = [-1,-1,-1,-1]
		self.readType = int(sysArgs[5])
		self.readFile_1 = sysArgs[1]
		self.readFile_2 = sysArgs[2]
		self.contigFile = sysArgs[3]
		self.plotPath = sysArgs[4]
		self.showPlot = int(sysArgs[6])
		
		
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
		regularPath = sys.argv[7]
		pathToBedTools = sys.argv[8]
		#regularPath = "/home/lexy/BWA_test/"
		core = CallerCore(sys.argv)
		core.doBWA(regularPath)
		core.doSAM(regularPath)
		core.doBAM(regularPath, pathToBedTools)
		core.doCoverageStats(regularPath)
		print core.getAllStats()
