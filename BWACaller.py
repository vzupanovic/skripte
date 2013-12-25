import subprocess
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import analyzer as stats
import operator
import matplotlib.ticker as mticker
from pylab import *


class BWACaller():
	def __init__(self, readType, referenceSequence, sequenceToMap1, sequenceToMap2):
		self.readType = readType
		self.referenceSequence = referenceSequence
		self.sequenceToMap1 = sequenceToMap1
		self.sequenceToMap2 = sequenceToMap2
		self.cmd1 = "bwa index "+str(self.referenceSequence)
		self.cmd2 = "bwa aln "+str(self.referenceSequence)+" "+str(self.sequenceToMap1)+" > sequence1.sai"
		self.cmd3 = "bwa aln "+str(self.referenceSequence)+" "+str(self.sequenceToMap2)+" > sequence2.sai"
		self.cmd4 = "bwa sampe "+str(self.referenceSequence)+" sequence1.sai sequence2.sai "+str(self.sequenceToMap1)+" "+str(self.sequenceToMap2)+" > alignment.sam"
		self.cmd5 = "bwa samse "+str(self.referenceSequence)+" sequence1.sai "+str(self.sequenceToMap1)+" > alignment.sam" 
	def calculateIndex(self):
		os.system(self.cmd1)
	def align(self):
		if (self.readType == 1):
			os.system(self.cmd2)
			os.system(self.cmd3)
		else:
			os.system(self.cmd2)
	def doSampe(self):
		os.system(self.cmd4)
	def doSamsa(self):
		os.system(self.cmd5)
		
class SAMTools():
	def __init__(self):
		self.cmd1 = "samtools view -bS alignment.sam | samtools sort - test_sorted"
		self.cmd2 = "samtools index test_sorted.bam test_sorted.bai"
		self.cmd3 = "samtools idxstats test_sorted.bam > analyze.txt"
		self.cmd4 = "samtools flagstat test_sorted.bam"
	def execute(self):
		os.system(self.cmd1)
		os.system(self.cmd2)
		os.system(self.cmd3)
		os.system(self.cmd4)
		
class BEDTools():
	def __init__(self, pathToBam, pathToCov): #dodaj da ima path!
		self.cmd1 = "~/bedtools/bin/bedtools genomecov -ibam "+pathToBam+"test_sorted.bam -dz > "+pathToCov+"test.bam.cov"
	def compute(self):
		os.system(self.cmd1)

class Coverage:
	def getContigData(self, inputData): #get basic contig info
		self.contiger = stats.ContigAnalyzer()
		self.contigs = []
		self.contigHeaders = []
		self.contigLengths = []
		self.contigs, self.contigHeaders, self.contigLengths = self.contiger.parseFasta(inputData)
		self.totalLenCont = self.contiger.getTotalLength()
		
	def getAlignmentData(self, inputData): #get alignment data
		stream = open(inputData, 'r')
		data = stream.readlines()
		alignData = {}
		self.uncoveredRegions = {}
		for header in self.contigHeaders:
			alignData[header] = []
			self.uncoveredRegions[header] = []
		#print alignData
		for line in data:
			line = line.strip()
			temp = []
			temp = line.split("\t")
			covValue = int(temp[-1])
			header = temp[0]
			alignData[header].append(covValue)
		self.alignData = alignData
		return alignData
		
	def getCoveragePerContig(self): #get all covered bases in contig, every base must be covered at least once
		coverageData = {}
		self.notCoveredContigs = []
		self.totalCoverage = 0
		print "[AN:] Getting coverage per contig..."
		for header in self.alignData:
			temp = self.alignData[header] 
			covered = 0
			for value in temp:
				if value > 0:
					covered += 1
			if len(temp) > 0:
				coverage = float(covered) / len(temp)
				coverageData[header] = coverage
				print "[AN:] Contig Id: ", header, " coverage: ", coverage
			else:
				self.notCoveredContigs.append(header)
				print "[AN:] Contig Id: ", header, " has 0 coverage"
		self.coverageData = coverageData
		print "[AN:] Contigs not covered by any read:"
		for contig in self.notCoveredContigs:
			print "\t",contig
		return coverageData
		
	def getUncovered(self):
		baseNum = 0
		allData = self.contiger.getContigDict()
		for header in self.notCoveredContigs:
			baseNum = baseNum + len(allData[header])
		return baseNum
			
		
	def getMaxCoveredContigs(self, n): #get contigs with highest number of bases covered
		sortedContigs = self.sortContigs()
		return sortedContigs[0:n]
			
	def sortContigs(self): #sort contigs by coverage
		sortedContigs = sorted(self.coverageData.iteritems(), key=operator.itemgetter(1), reverse = True)
		return sortedContigs
		
	def plotCoverage(self, contigId, path, showMe): #plot coverage of contig
		currentCovs = self.alignData[contigId]
		xAxis = np.arange(0, len(currentCovs))
		yAxis = np.array(currentCovs)
		plot(xAxis, yAxis, 'ro', markersize=3)
		xlabel('Relative position inside contig')
		ylabel('Number of reads')
		title('Contig coverage')
		grid(True)
		savefig(path+"/"+contigId+".png")
		if showMe == 1:
			show()
		
	def plotAllContigCov(self, path, show): #plot coverage of all contigs in .fasta file
		if not os.path.exists(path):
			os.makedirs(path)
		for header in self.alignData:
			if self.alignData[header] != []:
				print "Plotting:", header
				self.plotCoverage(header, path, show)
			else:
				print "Skipping contig:"+header+" none of the reads were mapped to it."				
		
	def getUncoveredRegions(self): #skip uncovered regions
		self.totalUncoveredBases = 0
		for header in self.alignData:
			temp = self.alignData[header] 
			for i in range(0, len(temp)):
				if temp[i] == 0:
					self.totalUncoveredBases += 1
					self.uncoveredRegions[header].append(i)
		print "Percentage of uncovered bases:",(float(self.totalUncoveredBases) + self.getUncovered()) / self.totalLenCont
		return self.totalUncoveredBases, float(self.totalUncoveredBases)/self.totalLenCont * 100
		
	def getRegionMaxCov(self, contigId): #get regions with highest coverage
		currentContigCov = self.alignData[contigId]
		if currentContigCov == []:
			return 0, None
		maxCoverage = max(currentContigCov)
		maxCovPos = [i for i, x in enumerate(currentContigCov) if x == maxCoverage]
		return maxCoverage, maxCovPos
		
	def getPotColapseRegions(self, n): #get potential colapse regions (high coverage)
		maxCovs = {}
		print "[AN:] Getting maximum coverage per contig..."
		for header in self.alignData:
			maxcov, pos = self.getRegionMaxCov(header)
			print "[AN:] Contig ID:",header, "max cov:", maxcov 
			maxCovs[header] = maxcov
		sortedCovs = sorted(maxCovs.iteritems(), key=operator.itemgetter(1), reverse = True)
		return sortedCovs[0:n]
		
class Initializer:
	def __init__(self, contigFile, readFile_1, readFile_2):
		self.contigFile = contigFile
		self.read1 = readFile_1
		self.read2 = readFile_2
