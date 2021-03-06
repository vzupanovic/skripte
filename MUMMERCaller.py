import subprocess
import os
import sys
import os.path
import re

class Mummer:
	def __init__(self, genomeFile, contigFile):
		self.genomeFile = genomeFile
		self.contigFile = contigFile
		self.refErrors = [-1,-1,-1,-1]
		self.qryErrors = [-1,-1,-1,-1]
		self.cmd1 = "nucmer -maxmatch -c 100 -p nucmer "+self.genomeFile+" "+self.contigFile
		self.cmd2 = "show-coords -r -c -l nucmer.delta > nucmer.coords"
		self.cmd3 = "show-snps -C nucmer.delta > nucmer.snps"
		self.cmd4 = "show-tiling nucmer.delta > nucmer.tiling"
		self.cmd5 = "delta-filter -m nucmer.delta > nucmer.delta.m"
		self.cmd6 = "mummerplot nucmer.delta.m"
		self.cmd7 = "dnadiff -d nucmer.delta"
		
	def doMummer(self):
		os.system(self.cmd1)
		os.system(self.cmd2)
		os.system(self.cmd3)
		os.system(self.cmd4)
		os.system(self.cmd5)
		
	def doPlot(self):
		os.system(self.cmd6)
		
	def getFeatures(self):
		os.system(self.cmd7)
		
	def printBasicMummerStats(self):
		stream = open('out.report', 'r')
		self.data = stream.readlines()
		for line in self.data:
			line = line.strip()
			print line
			
	def getErrors(self):
		for line in self.data:
			line = line.strip()
			line = "\t".join(line.split())
			if ("Breakpoints" in line):
				print '\t'+line
				temp = line.split("\t")
				self.refErrors[0] = temp[1]
				self.qryErrors[0] = temp[2]
			elif ("Relocations" in line):
				print '\t'+line
				temp = line.split("\t")
				self.refErrors[1] = temp[1]
				self.qryErrors[1] = temp[2]
			elif ("Translocations" in line):
				print '\t'+line
				temp = line.split("\t")
				self.refErrors[2] = temp[1]
				self.qryErrors[2] = temp[2]
			elif ("Inversions" in line):
				print '\t'+line
				temp = line.split("\t")
				self.refErrors[3] = temp[1]
				self.qryErrors[3] = temp[2]
				
	def getAllStats(self):
		return self.refErrors + self.qryErrors
			
		
if __name__ == "__main__":
	if len(sys.argv) < 3:
		print "[AN:] Usage: python MUMMERCaller.py genome_file.fasta contig_file.fasta"
	else:
		if os.path.isfile(sys.argv[1])!=True:
			print "[AN:] File "+str(sys.argv[1])+" doesn't exist. Exiting..."
			exit(-1)
		if os.path.isfile(sys.argv[2])!=True:
			print "[AN:] File "+str(sys.argv[2])+" doesn't exist. Exiting..."
			exit(-1)
		mummer = Mummer(sys.argv[1], sys.argv[2])
		print "[AN:] Running mummer..."
		mummer.doMummer()
		print "[AN:] Getting basic features..."
		mummer.getFeatures()
		print "[AN:] Printing basic stats..."
		mummer.printBasicMummerStats()
		print "[AN:] Getting the number of breakpoints, relocations, translocations and inversions (errors)[REF/QRY]...\n"
		mummer.getErrors()
		print "[AN:] Plotting..."
		mummer.doPlot()
		print mummer.getAllStats()
		
	
