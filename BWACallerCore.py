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
		self.readType = readType
		self.readFile_1 = readFile_1
		self.readFile_2 = readFile_2
		self.contigFile = contigFile
		self.plotPath = plotPath
		self.showPlot = showPlot
		
	def doBWA(self):
		self.bwa = BWACaller(self.readType, self.contigFile, self.readFile_1, self.readFile_2)
		self.bwa.calculateIndex()
		self.bwa.align()
		if self.readType == 1:
			self.bwa.doSamsa()
		else:
			self.bwa.doSampe()
		
	def doSAM(self):
		print "[AN:] Getting basic alignment stats..."
		self.samtools = SAMTools()
		self.samtools.execute()
		
	def doBAM(self):
		self.bedtools = BEDTools("~/celera_test/","~/celera_test/") #change this !
		self.bedtools.compute()
		
	def doCoverageStats(self):
		print "[AN:] Getting basic coverage stats..."
		coverage_stats = Coverage()
		coverage_stats.getContigData(self.contigFile)
		coverage_stats.getAlignmentData("test.bam.cov")
		coverage_stats.getCoveragePerContig()
		potCollaps = coverage_stats.getPotColapseRegions(2)
		print "[AN:] Potential colapse regions because of high coverage:"
		for contig in potCollaps:
			print "\tContig ID:",contig[0],"\tcoverage:",contig[1]
		print "[AN:] Plotting coverage per contig..."
		coverage_stats.plotAllContigCov(self.plotPath, self.showPlot)
		print "[AN:] Done..."
		

		

if __name__ == "__main__":
	if len(sys.argv) < 7:
		print "[AN:] Usage python BWACallerCore.py read1.fastq read2.fastq[null]\n\tcontig_file.fasta plot_folder read_type[0/1] show_plot[0/1]"
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
		core = CallerCore(readType, read1, read2, contigFile, plotFolder, showPlot)
		core.doBWA()
		core.doSAM()
		core.doBAM()
		core.doCoverageStats()
