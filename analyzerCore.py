import numpy as np
from matplotlib import pyplot as plt
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
import math
import sys
import os.path
from analyzer import *

class Core():
	def __init__(self, contigFile, scaffoldFile, qualitiyFile): #initialize analyzers
		self.contigFile = contigFile
		self.scaffoldFile = scaffoldFile
		self.qualityFile = qualityFile
		if contigFile != None:
			self.getContigStats()
			
		if scaffoldFile != None:
			self.getScaffoldStats()
			
		if scaffoldFile != None and contigFile != None:
			self.getMixedStats()
		
	def getContigStats(self): #get basic statistics for contigs
		self.contiger = ContigAnalyzer()
		print "\n[AN:] Getting basic statistics for contigs..."
		contigs, headers, lengths = self.contiger.parseFasta(self.contigFile)
		totalNum = self.contiger.getNumContigs()
		biggerThen = self.contiger.getCertainSize(25000)[0] #ulazni
		totalLen = self.contiger.getTotalLength()
		maxLen = self.contiger.getMaxContigLength()[1]
		minLen = self.contiger.getMinContigLength()[1]
		avgLen = self.contiger.getAvgContigLength()
		medLen = self.contiger.getMedContigLength() 
		eSize = self.contiger.getEsize(0) #ulazni param
		n25 = self.contiger.getNX(lengths, 25, self.contiger.getTotalLength())
		n50 = self.contiger.getNX(lengths, 50, self.contiger.getTotalLength())
		n56 = self.contiger.getNX(lengths, 56, self.contiger.getTotalLength())
		n75 = self.contiger.getNX(lengths, 75, self.contiger.getTotalLength())
		rc = self.contiger.getRatioRC(75) #ulazni param
		EM = self.contiger.getCorrelationEM(eSize, medLen)
		EN = self.contiger.getCorrelationEN(eSize, n50)
		print "[AN:] Total number of contigs: ", totalNum
		print "[AN:] Number of contigs bigger then size 25kB:", biggerThen
		print "[AN:] Total length of all contigs: ",totalLen
		print "[AN:] Maximum contig size: ", maxLen
		print "[AN:] Minimum contig size: ", minLen
		print "[AN:] Average contig length: ", avgLen
		print "[AN:] Median contig length: ",medLen
		print "[AN:] E-size with threshold 0: ", eSize #ulazni param!
		print "[AN:] N25 measure: ", n25
		print "[AN:] N50 measure: ", n50
		print "[AN:] N56 measure: ", n56
		print "[AN:] N75 measure: ", n75
		print "[AN:] Ratio between median and E-size: ", EM
		print "[AN:] Ratio between n50 and E-size: ", EN
		print "[AN:] Ratio between contig and read length: ", rc #ulazni param
		return [totalNum, biggerThen, totalLen, maxLen, minLen, avgLen, medLen, eSize, n25, n50, n56, n75]
		
	def getScaffoldStats(self): #get basic statistics for scaffolds
		self.scaffolder = ScaffoldAnalyzer()
		scalarStats = []
		print "\n[AN:] Getting basic statistics for scaffolds..."
		lengths, headers, scaffolds, totalSize = self.scaffolder.parseFasta(self.scaffoldFile)
		totalNum = self.scaffolder.getScaffNumber()
		maxSize = self.scaffolder.getMaxScaff()[1]
		minSize = self.scaffolder.getMinScaff()[1]
		avgSize = self.scaffolder.getAvgScaff()
		medSize = self.scaffolder.getMedScaffLength()
		certainNum = self.scaffolder.getCertainSizeNum(25000)[0] #ulazni
		n25 = self.scaffolder.getNX(lengths, 25, totalSize)
		n50 = self.scaffolder.getNX(lengths, 50, totalSize)
		n56 = self.scaffolder.getNX(lengths, 56, totalSize)
		n75 = self.scaffolder.getNX(lengths, 75, totalSize)
		eSize = self.scaffolder.getEsize(0) #ulazni
		EM = self.scaffolder.getCorrelationEM(eSize, medSize)
		EN = self.scaffolder.getCorrelationEN(eSize, n50)
		print "[AN:] Total number of scaffolds: ", totalNum
		print "[AN:] Maximum scaffold size: ", maxSize
		print "[AN:] Minimum scaffold size: ", minSize
		print "[AN:] Average scaffold size: ", avgSize
		print "[AN:] Median scaffold size: ", medSize
		print "[AN:] Total number of scaffolds bigger then size 25 kB", certainNum
		print "[AN:] N25 measure: ", n25
		print "[AN:] N50 measure: ", n50
		print "[AN:] N56 measure: ", n56
		print "[AN:] N75 measure: ", n75
		print "[AN:] Ratio between median and E-size: ", EM
		print "[AN:] Ratio between n50 and E-size: ", EN
		if (self.qualityFile == None):
			print "[AN:] Quality file not provided. Skipping..."
		else:
			print "[AN:] Average quality stats:"
			quals = self.scaffolder.parseQual(self.qualityFile)
			avgs, qual = self.scaffolder.getAverageQual(quals)
			for i in range(0, len(headers)):
				print "\t Scaffold: ", headers[i], " average quality: ",avgs[i]
			print "[AN:] Total average quality: ", qual
		return [totalNum, maxSize, minSize, avgSize, medSize, certainNum, n25, n50, n56, n75, EM, EN, avgs, qual]
		
	def getMixedStats(self): #get basic mixed statistics
		self.comparer = Comparer(self.contigFile, self.scaffoldFile)
		avgs = self.comparer.compAvgLength()
		n50s = self.comparer.compN50()
		nums = self.comparer.compNum()
		lens = self.comparer.compLenghts()
		maxs = self.comparer.compMaxLengths()
		print "\n[AN:] Getting basic mixed stats: "
		print "[AN:] Average scaffold length divied by average contig length: ",avgs
		print "[AN:] N50 of contigs divided by N50 of scaffolds: ", n50s
		print "[AN:] Number of contigs divided by number of scaffolds: ",nums
		print "[AN:] Total length of contigs divided by total length of scaffolds: ", lens
		print "[AN:] Maximal length of contig divided by maximal length of scaffold: ",maxs
		return [avgs, n50s, nums, lens, maxs]
				
if __name__ == "__main__":
	if len(sys.argv) < 2:
		print "[AN:] Usage: python analyzer.py path_to_contig_file path_to_scaffold_file[optional] quality_file[optional]"
	else:
		if len(sys.argv) == 2:
			contigFile = sys.argv[1]
			scaffoldFile = None
			qualityFile = None
			if os.path.isfile(sys.argv[1])!=True:
				print "[AN:] File "+str(sys.argv[1])+" doesn't exist. Exiting..."
				exit(-1)
		else:
			contigFile = sys.argv[1]
			scaffFile = sys.argv[2]
			qualityFile = None
			
			if os.path.isfile(sys.argv[1])!=True:
				print "[AN:] File "+str(sys.argv[1])+" doesn't exist. Exiting..."
				exit(-1)
			if os.path.isfile(sys.argv[2])!=True:
				print "[AN:] File "+str(sys.argv[2])+" doesn't exist. Exiting..."
				exit(-1)
				
			if (sys.argv) > 3:
				qualityFile = sys.argv[3]
				#print "DA"
				if os.path.isfile(sys.argv[3])!=True:
					print "[AN:] File "+str(sys.argv[3])+" doesn't exist. Exiting..."
					exit(-1)
				
		core = Core(contigFile, scaffFile, qualityFile)