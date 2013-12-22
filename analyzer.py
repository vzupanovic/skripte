import numpy as np
from matplotlib import pyplot as plt
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
import math

class ScaffoldAnalyzer:
	def parseFasta(self, fileName): #get scaffLength scaffIds scaffolds and totalSize
		self.scaffLenghts = []
		self.scaffIds = []
		self.scaffolds = []
		self.totalSize = 0
		for seqRecord in SeqIO.parse(fileName, "fasta"):
			self.scaffIds.append(seqRecord.id)
			self.scaffolds.append(seqRecord.seq)
			self.scaffLenghts.append(len(seqRecord))
			self.totalSize += len(seqRecord)
		print type(self.scaffolds[0])
		return self.scaffLenghts, self.scaffIds, self.scaffolds, self.totalSize
		
	def parseQual(self, fileName): #parse qual file, needed for SNP estimation
		qvals = []
		scaff = []
		stream = open(fileName, 'r')
		data = stream.readlines()
		for line in data:
			line = line.strip()
			if line[0] == '>' or line[0] == '@':
				scaff = ' '.join(scaff)
				qvals.append(scaff)
				scaff = []
			else:
				scaff.append(line)
		qvals = qvals[1:]
		scaff = ' '.join(scaff)
		qvals.append(scaff)
		return qvals
		
	def getAverageQual(self, qvals): #get average quality value
		splited = []
		avgs = []
		globalSum = 0
		globalLen = 0
		for block in qvals:
			splited = block.split(' ')
			s = 0
			for value in splited:
				s = s + int(value)
			avgs.append(float(s)/len(splited))
			globalSum += s
			globalLen += len(splited)
			
		return avgs, float(globalSum)/globalLen
				
		
	def mapPotSNP(self, inFileName, threshold, outFileName): #find potential SNPs (quality value)
		potSNPpos = {}
		SNPmatrix = {}
		SNPnum = 0
		blockNum = 0
		qvals = self.parseQual(inFileName)
		avgQual,globalQual = self.getAverageQual(qvals)
		writeStream = open(outFileName, 'w')
		for block in qvals:
			writeStream.write(self.scaffIds[blockNum]+"\n")
			SNProw = []
			blockSNP = []
			splited = block.split(' ')
			for i in range(0, len(splited)):
				if int(splited[i]) >= threshold:
					blockSNP.append(i)
					SNProw.append('X') #not row more block
					writeStream.write('X')
					SNPnum += 1
				else:
					SNProw.append('O')
					writeStream.write('O')
			potSNPpos[self.scaffIds[blockNum]] = blockSNP
			SNPmatrix[self.scaffIds[blockNum]] = SNProw
			blockNum = blockNum + 1
			writeStream.write('\n')
		writeStream.close()
		return potSNPpos, SNPmatrix, SNPnum, avgQual, globalQual
			
				
	def getEstimatedSize(self):
		estimatedSize = 0 #estimated genome size is equal to number of bases in scaffold
		for scaffold in self.scaffolds:
			estimatedSize = estimatedSize + len(scaffold)
		return estimatedSize
		
	def getScaffNumber(self): #get number of scaffolds
		return len(self.scaffolds)
		
	def sortScaff(self): #sort scaffolds by length
		new = sorted(self.scaffolds, key = len)
		self.scaffolds = new[:]
		
	def getContigCoverage(self, contigsFile):
		contigsAnalyze = ContigAnalyzer()
		contigs, headers, contigLenghts = contigsAnalyze.parseFasta(contigsFile)
		contigTotalLen = contigsAnalyze.getTotalLength(contigs)
		coverage = float(contigTotalLen) / self.totalSize * 100
		return coverage
		
	def getMaxScaff(self): #maximum scaffold size
		largestScaff = max(self.scaffolds, key=len)
		return largestScaff, len(largestScaff)
		
	def getMinScaff(self): #minimum scaffold size
		largestScaff = min(self.scaffolds, key=len)
		return largestScaff, len(largestScaff)
		
	def getAvgScaff(self): #average scaffold size
		totalNum = len(self.scaffolds)
		lengths = []
		totalLen = 0
		for scaffold in self.scaffolds:
			totalLen = totalLen + len(scaffold)
			lengths.append(len(scaffold))
		avg = totalLen/totalNum
		return avg
		
	def getCertainSizeNum(self, size): #get number of scaffolds bigger then certain size
		goodScaff = []
		for scaffold in self.scaffolds:
			if len(scaffold) >= size:
				goodScaff.append(scaffold)
		return len(goodScaff), goodScaff
		
	def getNX(self, lenghts, measure, totalSize): #n25 n50 n75 nX calculation
		orderdScaffLength = sorted(lenghts, reverse = True)
		NX = float(measure) / 100 * totalSize
		sumLenght = 0
		measureValue = -1
		for lenght in orderdScaffLength:
			sumLenght += lenght
			if sumLenght >= NX:
				measureValue = lenght
				break
		return measureValue
		
	def getEsize(self, threshold): #eSize calculation
		genomeSize = self.getEstimatedSize()
		Esize = 0
		for length in self.scaffLenghts:
			if (length > threshold):
				Esize += math.pow(length, 2) / genomeSize
		return Esize
		
	def getCorrelationEN(self, Esize, N50): #correlation between n50 and eSize
		return float(Esize/N50)	
		
	def getLenList(self): #get sizes of all scaffolds
		return self.scaffLenghts
		
	def getMedScaffLength(self): #median calculation
		self.sortScaff()
		scaffNum = len(self.scaffolds)
		if (scaffNum % 2 == 1):
			median = len(self.scaffolds[(scaffNum + 1) / 2])
		else:
			median = (len(self.scaffolds[(scaffNum)/ 2 - 1]) + len(self.scaffolds[(scaffNum) / 2])) / 2
		return median
		
	def getCorrelationEM(self, eSize, median): #corelation between median and eSize
		return float(eSize)/median
			
		
class ContigAnalyzer:
	def parseFasta(self, fileName): #parse input file to get contigs, contigs id-s and contig lenghts
		self.contigLengths = []
		self.contigIds = []
		self.contigs = []
		for seqRecord in SeqIO.parse(fileName, "fasta"):
			self.contigIds.append(seqRecord.id)
			self.contigs.append(seqRecord.seq)
			self.contigLengths.append(len(seqRecord))
		return self.contigs, self.contigIds, self.contigLengths
		
	def getNumContigs(self): #get number of contigs
		return len(self.contigLengths)	
			
	def getAvgContigLength(self): #get average contig length
		totalNum = len(self.contigs)
		lengths = []
		totalLen = 0
		for contig in self.contigs:
			totalLen = totalLen + len(contig)
		avg = float(totalLen)/totalNum
		return avg
		
	def getRatioRC(self, readLength): #get ratio between read and contig length
		avgLength = self.getAvgContigLength()
		return float(avgLength)/readLength
		
	def getAvgMedDiff(self): #get difference between median and average contig length
		avgLength = self.getAvgContigLength()
		medLength = self.getMedContigLength()
		return abs(avgLength - medLength)
		
	def orderContigs(self): #sort contigs by length
		new = sorted(self.contigs, key = len)
		self.contigs = new[:]
		
	def getLenList(self): #get list of contigs length
		return self.contigLengths
		
	def getMedContigLength(self): #get median contig length
		self.orderContigs()
		contigsNum = self.getNumContigs()
		if (contigsNum % 2 == 1):
			median = len(self.contigs[(contigsNum + 1) / 2])
		else:
			median = (len(self.contigs[(contigsNum)/ 2 - 1]) + len(self.contigs[(contigsNum) / 2])) / 2
		return median
			
	def getMaxContigLength(self): #get biggest contig and its length
		largestContig = max(self.contigs, key=len)
		return largestContig, len(largestContig)
		
	def getMinContigLength(self): #get smallest contig and its length
		smalestContig = min(self.contigs, key=len)
		return smalestContig, len(smalestContig)
		
	def getTotalLength(self): #get sum of contig lengths
		totalLen = 0
		for contig in self.contigs:
			totalLen = totalLen + len(contig)
		return totalLen
		
	def getCertainSize(self, size): #get number of contigs bigger then certain size
		goodContigs = []
		for contig in self.contigs:
			if len(contig) >= size:
				goodContigs.append(contig)
		return len(goodContigs), goodContigs
		
	def getEsize(self,threshold): #get eSize
		genomeSize = self.getTotalLength()
		Esize = 0
		for length in self.contigLengths:
			if (length > threshold):
				Esize += math.pow(length, 2) / genomeSize
		return Esize
		
	def getNX(self, lenghts, measure, totalSize): #n25 n50 n75 nX
		orderdContigLength = sorted(lenghts, reverse = True)
		NX = float(measure) / 100 * totalSize
		sumLenght = 0
		measureValue = -1
		for lenght in orderdContigLength:
			sumLenght += lenght
			if sumLenght >= NX:
				measureValue = lenght
				break
		return measureValue
		

class Comparer():
	def __init__(self, inputFileContigs, inputFileScaff):
		self.scaffolder = ScaffoldAnalyzer()
		self.contiger = ContigAnalyzer()
		self.scaffolder.parseFasta(inputFileContigs)
		self.contiger.parseFasta(inputFileScaff)
		
		
	def compAvgLength(self): #compare average length of contig to average length of scaffold
		contigLen = self.contiger.getAvgContigLength()
		scaffLen = self.scaffolder.getAvgScaff()
		return float(scaffLen) / contigLen
		
	def compN50(self): #compare n50 measures
		totalSizeCont = self.contiger.getTotalLength()
		totalSizeScaff = self.scaffolder.getEstimatedSize()
		sizeListCon = self.contiger.getLenList()
		sizeListScaff = self.scaffolder.getLenList()
		n50s = self.scaffolder.getNX(sizeListScaff, 50, totalSizeScaff)
		n50c = self.contiger.getNX(sizeListCon, 50, totalSizeCont)
		return float(n50c) / n50s
		
	def compNum(self): #compare total number of contigs to total number of scaffolds
		scaffNum = self.scaffolder.getScaffNumber()
		contigNum = self.contiger.getNumContigs()
		return float(contigNum)/scaffNum
		
	def compLenghts(self): #compare total length of contig file to total length of scaffold file
		totalLenCont = self.contiger.getTotalLength()
		totalLenScaff = self.scaffolder.getEstimatedSize()
		return float(totalLenCont) / totalLenScaff
	
	def compMaxLengths(self): #compare maximal contig to maximal scaffold
		largestCont, maxLenCont = self.contiger.getMaxContigLength()
		largestScaff, maxLenScaff = self.scaffolder.getMaxScaff()
		return float(maxLenCont) / maxLenScaff
		
		
#if __name__ == "__main__":
#	scaffolder = ScaffoldAnalyzer()
#	contiger = ContigAnalyzer()
#	comparer = Comparer("contigzi.fasta", "scaffolds.fasta")
	#do something...

