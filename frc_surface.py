################################################################
# Svaki FRC file 3 kolone kolona 1: feature threshold T,
# kolona 2: n50 contiga u zavisnosti od T iz kolone 1
# kolona 3: kumulativna velicina contiga ciji je broj
# znacajki <= T
# (thrs, n50, contig_coverage * 100.00) trojke (FRC.txt)

import numpy as np
from matplotlib import pyplot as plt
import os

class Curve:
	def approximateCurve(self, inputDat):
		featureThrs = []
		appCov = []
		n50s = []
		featureFile = open("FRC.txt",'r')
		features = featureFile.readlines()
		for feature in features:
			feature = feature.strip()
			data = []
			data = feature.split()
			print data
			featureThrs.append(float(data[0]))
			n50s.append(float(data[1]))
			appCov.append(float(data[2]))
		featureThrs = np.array(featureThrs)
		appCov = np.array(appCov)
		n50s = np.array(n50s)
		coef = np.polyfit(featureThrs, appCov, 3)
		polyObj = np.poly1d(coef)
		print polyObj, coef
		return coef, polyObj, appCov, featureThrs
		
	def plotCurve(self, inputData, coef, p, appCov, featureThrs, path): #p je polyObj
		if not os.path.exists(path):
			os.makedirs(path)
		xp = np.linspace(0, 2000)
		print appCov, featureThrs
		plt.ylim(0, 120)
		plt.plot(featureThrs, appCov, xp, p(xp))
		plt.xlabel('Relative position inside contig')
		plt.ylabel('Number of reads')
		plt.title('Contig coverage')
		plt.grid(True)
		plt.savefig(path+"frc.png")
		plt.show()
		
	def getSurface(self, appCov, featureThrs):
		return np.trapz(featureThrs, appCov)
		
if __name__ == "__main__":
	curve = Curve()
	coef, polyObj, appCov, featureThrs = curve.approximateCurve("FRC.txt")
	curve.plotCurve("FRC.txt", coef, polyObj, appCov, featureThrs, "nesto/")
	print curve.getSurface(appCov, featureThrs)
	
