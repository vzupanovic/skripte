import numpy as np
from matplotlib import pyplot as plt
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
import math
import sys
import os.path

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
		
	def getCorrelationEN(self, Esize, N50): #correlation between n50 and eSize
		return float(Esize/N50)	
		
	def getCorrelationEM(self, eSize, median): #corelation between median and eSize
		return float(eSize)/median
		

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
		
class Core():
	def __init__(self, contigFile, scaffoldFile, qualitiyFile):
		self.contigFile = contigFile
		self.scaffoldFile = scaffoldFile
		self.qualityFile = qualityFile
		if contigFile != None:
			self.getContigStats()
			
		if scaffoldFile != None:
			self.getScaffoldStats()
			
		if scaffoldFile != None and contigFile != None:
			self.getMixedStats()
		
	def getContigStats(self):
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
		
	def getScaffoldStats(self):
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
		
		
	def getMixedStats(self):
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


