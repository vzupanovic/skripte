import subprocess
import os
import numpy as np
import math
import matplotlib.pyplot as pltS
class BWACaller():
	def __init__(self, readType, referenceSequence, sequenceToMap1, sequenceToMap2):
		self.readType = readType
		self.referenceSequence = referenceSequence
		self.sequenceToMap1 = sequenceToMap1
		self.sequenceToMap2 = sequenceToMap2
		self.cmd1 = "bwa index "+str(self.referenceSequence)
		self.cmd2 = "bwa aln "+str(self.referenceSequence)+" "+str(self.sequenceToMap1)+" > sequence1.sai"
		self.cmd3 = "bwa aln "+str(self.referenceSequence)+" "+str(self.sequenceToMap2)+" > sequence2.sai"
		self.cmd4 = "bwa sampe "+str(self.referenceSequence)+" sequence1.sai sequence2.sai "+str(self.sequenceToMap1)+" "+str(self.sequenceToMap2)+" > out.sam"
		self.cmd5 = "bwa samse "+str(self.referenceSequence)+" sequence1.sai "+str(self.sequenceToMap1)+" > alignment.sam"
		#self.cmd6 = "bwa 
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
	def __init__(self):
		self.cmd1 = "~/bedtools/bin/bedtools genomecov -ibam ~/probaj2/test_sorted.bam -d > ~/probaj2/guzica.bam.cov"
	def compute(self):
		os.system(self.cmd1)

class Coverage:
	def parseInput(self, inputData):
		stream = open(inputData, 'r')
		data = stream.readlines()
		counter = 0
		self.xvalues = np.arange(1, len(data) + 1)
		self.yvalues = []
		#print type(self.xvalues)
		#print self.xvalues
		for line in data:
			temp = []
			line = line.strip()
			temp = line.split("\t")
			self.yvalues.append(int(temp[-1]))
			if int(temp[-1]) > 0:
				counter = counter + 1
		self.yvalues = np.array(self.yvalues)
		#print self.yvalues, type(self.yvalues)
		self.dataLen = len(data)
		self.counter = counter
		if (self.dataLen > 0):
			self.totalAverageCoverage = float(counter) / self.dataLen #each position covered at least once
		else:
			self.totalAverageCoverage = 0
		#print len(data), counter, float(counter)/len(data)
		
	def outputAverageCoverage(self):
		print "Genome coverage estimated to: ",self.totalAverageCoverage*100,"%"
		
	def outputCoveragePerPosition(self):
		print "Position:\tCoverage:"
		for i in range(0, len(self.yvalues)):
			print self.xvalues[i], " => \t", self.yvalues[i]
	
	def plotHistogram(self):
		plt.bar(self.xvalues, self.yvalues, color = 'r')
		#plt.bar((0,1,2,3,4), (0,1,4,2,0))
		plt.show()	
	
	
if __name__ == "__main__":
	bwa = BWACaller(0, "referenca.fasta", "contigi.fasta", "null")
	bwa.calculateIndex()
	bwa.align()
	bwa.doSamsa()
	samtools = SAMTools()
	samtools.execute()
	bedtools = BEDTools()
	bedtools.compute()
	coverage_stats = Coverage()
	coverage_stats.parseInput("guzica.bam.cov")
	coverage_stats.outputAverageCoverage()


			
		
	
		
		

