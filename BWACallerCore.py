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
	def __init__(self, readFile1, readFile2, contigFile, plotPath, readType, showPlot, regularPath, pathToBedTools):
		self.maxCovReg = [-1,-1,-1,-1]
		
		self.readFile_1 = readFile1;			# Path to the first reads file.
		self.readFile_2 = readFile2;			# Path to the second reads file.
		self.contigFile = contigFile;			# Path to the contig file that resulted from the above two files.
		self.plotPath = plotPath;			# Path for output plots.
		self.readType = int(readType);			# readType='0' for single end, and '1' for paired end reads.
		self.showPlot = int(showPlot);			# showPlot='0' if plots are not supposed to be displayed right away, otherwise showplot='1'.
		self.regularPath = regularPath;			# regularPath is the path for output alignment files.
		self.pathToBedTools = pathToBedTools;		# Path to BED tools.
		
		if os.path.isfile(self.readFile_1)!=True:
			print "[AN:] File "+str(self.readFile_1)+" doesn't exist. Exiting..."
			exit(-1)
		if os.path.isfile(self.readFile_2)!=True:
			print "[AN:] File "+str(self.readFile_2)+" doesn't exist. Exiting..."
			exit(-1)
		if os.path.isfile(self.contigFile)!=True:
			print "[AN:] File "+str(self.contigFile)+" doesn't exist. Exiting..."
			exit(-1)
		
	def doBWA(self):
		self.bwa = BWACaller(self.readType, self.contigFile, self.readFile_1, self.readFile_2, self.regularPath)
		self.bwa.calculateIndex()
		self.bwa.align()
		if self.readType == 1:
			self.bwa.doSamsa()
		else:
			self.bwa.doSampe()
		
	def doSAM(self):
		print "[AN:] Getting basic alignment stats..."
		self.samtools = SAMTools(self.regularPath)
		self.samtools.execute()
		
	def doBAM(self):
		self.bedtools = BEDTools(self.regularPath, self.pathToBedTools) #change this !
		self.bedtools.compute()
		
	def doCoverageStats(self):
		print "[AN:] Getting basic coverage stats..."
		coverage_stats = Coverage()
		coverage_stats.getContigData(self.contigFile)
		coverage_stats.getAlignmentData(self.regularPath+"test.bam.cov")
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
#		regularPath = sys.argv[7]
#		pathToBedTools = sys.argv[8]
		#regularPath = "/home/lexy/BWA_test/"
		
		readFile_1 = sys.argv[1]
		readFile_2 = sys.argv[2]
		contigFile = sys.argv[3]
		plotPath = sys.argv[4]
		readType = int(sys.argv[5])
		showPlot = int(sys.argv[6])
		regularPath = sys.argv[7]
		pathToBedTools = sys.argv[8]
		
		core = CallerCore(readFile_1, readFile_2, contigFile, plotPath, readType, showPlot, regularPath, pathToBedTools)
		core.doBWA()
		core.doSAM()
		core.doBAM()
		core.doCoverageStats()
		print core.getAllStats()
