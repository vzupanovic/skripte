#############################################################
# Basic viewer version 1.0
#: Simple tool used for the data visualisation 
#: Plots various assembly statistics
#: Input: .csv and .p files
#: Usage: python benchmarkGUI.py
#############################################################

import wx
import gui
from math import *
import numpy as np
import os
from threading import Thread
import matplotlib
matplotlib.use('WXAgg')
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
import summaryParser
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
		FigureCanvasWxAgg as FigCanvas, \
		NavigationToolbar2WxAgg as NavigationToolbar
import pickle
from matplotlib import rcParams
import subprocess
from matplotlib.ticker import ScalarFormatter
import time


class BenchmarkMainFrame(gui.MainFrame):
	def __init__(self, parent):
		gui.MainFrame.__init__(self, parent)
		self.summaryFiles = []
		self.activeSummaryFiles = []
		self.testedAssemblers = []
		self.summaryParsers = []
		self.summaryLabels = []
		self.covData = {}
		self.covDataKeys = []
		self.covDataValues = []
		self.plotIndex = 0
		self.plotDir = ''
		self.newEntry = ''
		self.xAttribute = 'l'
		self.yAttribute = 'cN50'
		self.xUnits = 'bp'
		self.yUnits = 'bp'
		self.xScale = 'linear'
		self.yScale = 'linear'
		self.readFile = ''
		self.referencePickerName = ''
		self.contigPickerName = ''
		self.deBrujinAssemblers = ['abyss','ray','soap','velvet']
		self.deBrujin = ('ABySS assembler', 'SOAPdenovo2', 'Velvet', 'Ray assembler')
		self.styles = []
		rcParams.update({'figure.autolayout': True})
		self.atributes = ['l', 'cov', 'N', 'd', 'e', 'r', 'R', 'X', 'A', 'D']
		self.units = ['bp', 'coverage', 'num reads', '', '', '', '', '', '', '']
		self.detailsDict = {'cTotalNum':'(number of contigs)', 'cBiggerThen':'(num. of contigs bigger then s)',
							'cTotalLen' : '(total length of contigs)', 'cMaxLen' : '(maximum length of contigs)',
							'cMinLen' : '(minimum length of contigs)', 'cAvgLen' : '(average length of contigs)',
							'cMedLen' : '(median length of contigs)', 'cN50':'(N50 size of contigs)',
							'cN25':'(N25 size of contigs)', 'cN75':'(N75 size of contigs)', 'cN56':'(N56 size of contigs)',
							'sTotalNum':'(number of scaffolds)', 'sBiggerThen':'(num. of scaffolds bigger then s)',
							'sTotalLen' : '(total length of scaffolds)','sMaxSize' : '(maximum length of scaff.)',
							'sMinSize' : '(minimum length of scaff.','sAvgLen' : '(average length of scaff.)',
							'sMedSize' : '(median length of scaff.)', 'sN50':'(N50 size of scaffolds)',
							'sN25':'(N25 size of scaffolds)','sN56':'(N56 size of scaffolds)',
							'sN75':'(N75 size of scaffolds)','sEM':'(ratio between median and E-size[scaff.])',
							'sEN':'(ratio between n50 and E-size)','mAvgs':'(average length of scaff./average length of cont.)',
							'mN50s':'(N50[contigs]/N50[scaffolds])','mNums':'([number of contigs]/[number of scaffolds])',
							'mLens':'([total len. of cont.]/[total len. of scaff.])', 'mMaxs':'([max length of con.]/[max length of scaff.])',
							'totalRealTime':'(total execution time of all steps of the assembly process)',
							'totalCpuTime':'(total CPU time of all steps of the assembly process)',
							'totalRSS':'(peak memory usage [Resident Set Size])',
							'l':'(read length)',
							}
		self.atributes += ['totalRealTime', 'totalCpuTime', 'totalRSS', 'totalPSS', 'totalVmSize', 'cTotalNum', 'cBiggerThen', \
						   'cTotalLen', 'cMaxLen', 'cMinLen', 'cAvgLen', 'cMedLen', 'cESize', 'cN25', 'cN50', 'cN56', 'cN75', 'sTotalNum', \
						   'sMaxSize', 'sMinSize', 'sAvgSize', 'sMedSize', 'sCertainNum', 'sN25', 'sN50', 'sN56', 'sN75', 'sEM', 'sEN', 'sQual', 'mAvgs', 'mN50s', 'mNums', 'mLens', 'mMaxs','cReferenceCoverage']
		self.units += ['sec', 'sec', 'MB', 'MB', 'MB', '', '', \
					   'bp', 'bp', 'bp', 'bp', 'bp', '', 'bp', 'bp', 'bp', 'bp', '', \
					   'bp', 'bp', 'bp', 'bp', '', 'bp', 'bp', 'bp', 'bp', '', '', \
					   '', '', '', '', '', '',''];
		self.selectionMap = {0:'linear',1:'logarithmic'}
		self.checkListBoxInitialItems = {"abyss":0,"minimus":1,"sga":2,"pasqual":3,"soap":4,"celera":5,"velvet":6,"readjoiner":7, "ray":8}
		self.assemblerTypes = {}
		print "[AN:] Basic viewer 1.0 started at: ", time.strftime("%H:%M:%S")
		self.CreatePlot()
		self.createReadsPlot()
		
	def handleCheckListBox(self):
		checkedIndeces = self.choiceList.GetChecked()
		checkedByName = self.choiceList.GetCheckedStrings()
		self.styles = []
		if self.choiceBrujin.GetValue() and self.newEntry != '':
			self.deBrujinAssemblers.append(self.newEntry.split('-')[-3])
		
		self.summaryLabels = checkedByName[:]
		self.activeSummaryFiles = []
		for ind in checkedIndeces:
			selectedString = self.checkListBoxInitialItems.keys()[self.checkListBoxInitialItems.values().index(ind)]
			if self.checkAssemblerType(selectedString):
				self.styles.append('--')
			else:
				self.styles.append('-')
			for fileName in self.summaryFiles:
				if fileName.find(selectedString) != -1:
					self.activeSummaryFiles.append(fileName)
				
		
		
	def addNewEntryToListBox(self, event):
		if self.newEntry != '':
			newAssembler = self.newEntry.split('-')[-3]
			self.choiceList.Append(newAssembler)
			self.checkListBoxInitialItems[newAssembler] = len(self.checkListBoxInitialItems)
			print self.checkListBoxInitialItems[newAssembler]
		else:
			dial = wx.MessageDialog(None, "Path to the new assembly data is not specified!\nPlease specify the path to the missing summary file.", 'Adding new assembly data error',  wx.ICON_ERROR)
			dial.ShowModal()


	def generatePlot(self):
		print "[AN:] Generating plot..."
		self.axes.clear()
		title = "[" + self.xAttribute + "/" + self.yAttribute + "]" + " dependency:" 
		self.xUnits = self.units[self.atributes.index(self.xAttribute)]
		self.yUnits = self.units[self.atributes.index(self.yAttribute)]
		valuesX = []
		valuesY = []
		allY = []
		counter = 0
		for summary in self.summaryParsers:
			plotMe = False
			valuesX = summary.getValuesForParam(self.xAttribute)
			valuesY = summary.getValuesForParam(self.yAttribute)
			allY += valuesY
			valuesY = np.array(valuesY).astype(np.double)
			valuesYPoints = valuesY

			goodX = []
			goodY = []
			missingX = []
			missingXIndex = []

			for m in range(0, len(valuesY)):
				if valuesY[m] != -1.0:
					plotMe = True
					goodX.append(valuesX[m])
					goodY.append(valuesY[m])
				else:
					missingX.append(valuesX[m])
					missingXIndex.append(m)
			missingY = []
			if len(goodX)> 0:
				missingY = np.interp(missingX, goodX, goodY)
				plotMe = True

			if plotMe == True:
				for m in range(0, len(missingXIndex)):
					valuesY[missingXIndex[m]] = missingY[m]
			
			currentColor = cm.jet(1.*counter/len(self.summaryParsers))
			self.axes.plot(valuesY, ls = self.styles[counter], label=self.summaryLabels[counter], color=currentColor)
			self.axes.plot(valuesYPoints, 'o', color=currentColor,  markersize=4)
			counter = counter + 1
	   
		self.setScales()
		
		handles, labels = self.axes.get_legend_handles_labels()
		leg = self.axes.legend(handles[::-1], labels[::-1], loc='best',prop={'size':10},fancybox=True)
		leg.get_frame().set_alpha(0.5)
		self.axes.axhline(0, color='k')
		
		self.axes.set_xticks(range(0, len(valuesX)), minor=False)
		self.axes.set_xticklabels(valuesX, rotation=70)

		if self.xAttribute in self.detailsDict:
			infoX = self.detailsDict[self.xAttribute]
		else:
			infoX = ''
		if self.yAttribute in self.detailsDict:
			infoY = self.detailsDict[self.yAttribute]
		else:
			infoY = ''
		self.axes.set_xlabel(self.xAttribute + infoX + ' ' + self.xUnits,fontsize=10)
		self.axes.set_ylabel(self.yAttribute + infoY + ' ' + self.yUnits,fontsize=10)
		allY = sorted(list(set([float(value) for value in allY])))
		self.axes.yaxis.set_major_locator(MaxNLocator(nbins=15, prune='upper'))
		self.axes.yaxis.set_minor_locator(MaxNLocator(nbins=15, prune='upper'))
		
		self.axes.grid()
		self.canvas.draw()
		
		print "[AN:] Done..."
	

	def updateBasicStatsView(self):
		if self.activeSummaryFiles:
			try:
				self.generatePlot()
			except Exception:
				dial = wx.MessageDialog(None, "Something went wrong, most probably you are trying to use log scale on negative values, check your input and try to plot again.", 'Plotting error', wx.ICON_ERROR)
				dial.ShowModal()
				
		else:
			dial = wx.MessageDialog(None, "Please, load some data to plot, there is nothing to plot at the moment.", 'Info', wx.OK)
			dial.ShowModal()


	def setScales(self):
		if self.xScale == 'linear':
			self.axes.set_xscale('linear')
		if self.xScale == 'logarithmic':
			self.axes.set_xscale('log')
		if self.yScale == 'linear':
			self.axes.set_yscale('linear')
		if self.yScale == 'logarithmic':
			self.axes.set_yscale('log')
		

	def CreatePlot(self): #just sample plot, it will be replaced with real one after you load some data
		formatter = ScalarFormatter()
		formatter.set_scientific(True)
		formatter.set_powerlimits((0,0)) 
		self.figure = Figure()
		self.figure.set_facecolor('white')
		self.figure.subplots_adjust(bottom=0.3, left=0.25) 
		self.axes = self.figure.add_subplot(111)
		self.axes.xaxis.set_major_formatter(formatter) 
		self.axes.yaxis.set_major_formatter(formatter) 
		x = np.arange(0,6,.01)
		y = np.sin(x**2)*np.exp(-x)
		self.axes.plot(x,y, ls = 'dotted',label = "This is just a sample plot and it will be replaced with\nthe real plot once when you load some data...")
		self.setScales()
	
		handles, labels = self.axes.get_legend_handles_labels()
		self.axes.legend(handles[::-1], labels[::-1], fancybox=True)
		frame=self.axes.get_frame()
		frame.set_alpha(0.4) 
		self.canvas = FigCanvas(self.plotPanel, wx.ID_ANY, self.figure) #jako bitna stavka
		return 1

	def createReadsPlot(self):
		currentCovs = [4,5,6,7,4,5,2]
		self.readFigure = Figure(figsize=(3.2, 2.2), dpi=150)
		self.readFigure.set_facecolor('white')
		self.readAxes = self.readFigure.add_subplot(111)
		self.readAxes.xaxis.set_tick_params(labelsize=5)
		self.readAxes.yaxis.set_tick_params(labelsize=5)
		x = np.arange(0, len(currentCovs))
		y = np.array(currentCovs)
		self.readAxes.plot(x,y, label = "test")
		self.readFigure.subplots_adjust(bottom=0.25, left = 0.15) 
		
		self.readCanvas = FigCanvas(self.readPanel, wx.ID_ANY, self.readFigure)

	def updateReadsPlotNext(self, event):
		print "[AN:] Plotting statistics for next contig..."
		self.plotIndex = self.plotIndex + 1
		if (self.plotIndex >= len(self.covDataKeys)):
			self.plotIndex = self.plotIndex - 1
			dial = wx.MessageDialog(None, "Nothing more to plot...", 'Info', wx.OK)
			dial.ShowModal()
			print "[AN:] Aborted..."
		else:
			self.updateReadView()
			print "[AN:] Done..."

	def updateReadsPlotPrevious(self, event):
		print "[AN:] Plotting statistics for previous contig..."
		self.plotIndex = self.plotIndex - 1
		if (self.plotIndex < 0):
			self.plotIndex = self.plotIndex + 1
			dial = wx.MessageDialog(None, "Can't go back, you are already on the first plot...", 'Info', wx.OK)
			dial.ShowModal()
			print "[AN:] Aborted..."
		else:
			self.updateReadView()
			print "[AN:] Done..."
	

	def getReadsData(self, event):
		print "[AN:] Loading mapping data..."
		self.plotIndex = 0
		self.readFile = self.readsDataPicker.GetPath()
		self.parseReadFile()
		self.readPlotName.SetLabel(self.covDataKeys[self.plotIndex])
		print "[AN:] Done..."
		self.updateReadView()

	def parseReadFile(self):
		self.covData = pickle.load(open(self.readFile, "rb"))
		self.covDataKeys = list(self.covData.keys())
		self.covDataValues = list(self.covData.values())

	def updateReadView(self):
		self.readPlotName.SetLabel(self.covDataKeys[self.plotIndex])
		self.readAxes.clear()
		currentCovs = self.covDataValues[self.plotIndex]
		x = np.arange(0, len(currentCovs))
		y = np.array(currentCovs)
		self.readAxes.plot(x,y, label = "test")
		self.readAxes.set_xlabel(self.covDataKeys[self.plotIndex]+"-(location on the contig)",fontsize=5)
		self.readAxes.set_ylabel('numReads-(number of mapped reads)',fontsize=5)
		self.readAxes.grid()
		self.readCanvas.draw()
		
	def checkAssemblerType(self, name):
		for assembler in self.deBrujinAssemblers:
			if name.find(assembler) != -1:
				return True
		return False
		
	def selectNewSummaryFile(self, event):
		self.newEntry = self.newAssemblerPicker.GetPath()
		self.summaryFiles.append(self.newEntry)


	def setSummaryDirectory(self, event):
		self.summaryLabels = []
		self.summaryDir = self.summaryDirPic.GetPath()
		print "[AN:] Loading summary directory: ",self.summaryDir
		self.summaryFiles = os.listdir(self.summaryDir)[:]
		temp = []
		for f in self.summaryFiles:
			if f.find("summary") != -1 and f.endswith('.csv'):
				print "\t",f
				self.summaryLabels.append(f)
				if self.checkAssemblerType(f):
					self.styles.append('--')
				else:
					self.styles.append('-')
			f = os.path.join(self.summaryDir, f)
			temp.append(f)
		self.summaryFiles = temp[:]
		validFiles = 0
		validFilesList = []
		for summary in self.summaryFiles:
			if summary.find("summary") != -1 and summary.endswith('.csv'):
				validFiles += 1
				validFilesList.append(summary)
		if validFiles == 0:
			dial = wx.MessageDialog(None, "Selected directory doesn't contain any summary files!", 'Info', wx.OK)
			dial.ShowModal()
		self.summaryFiles = validFilesList[:] #important, list of summary files!
		self.setInitialCheckBoxes() #after loading dir setting initial checkboxes
		print "[AN:] Done..."

	def setInitialCheckBoxes(self):
		selectedIndeces = []
		for summary in self.summaryFiles:
			if summary.find("abyss") != -1:
				selectedIndeces.append(self.checkListBoxInitialItems["abyss"])
			if summary.find("celera") != -1:
				selectedIndeces.append(self.checkListBoxInitialItems["celera"])
			if summary.find("minimus") != -1:
				selectedIndeces.append(self.checkListBoxInitialItems["minimus"])
			if summary.find("velvet") != -1:
				selectedIndeces.append(self.checkListBoxInitialItems["velvet"])
			if summary.find("sga") != -1:
				selectedIndeces.append(self.checkListBoxInitialItems["sga"])
			if summary.find("readjoiner") != -1:
				selectedIndeces.append(self.checkListBoxInitialItems["readjoiner"])
			if summary.find("pasqual") != -1:
				selectedIndeces.append(self.checkListBoxInitialItems["pasqual"])
			if summary.find("ray") != -1:
				selectedIndeces.append(self.checkListBoxInitialItems["ray"])
			if summary.find("soap") != -1:
				selectedIndeces.append(self.checkListBoxInitialItems["soap"])
		self.choiceList.SetChecked(selectedIndeces)


	def clearAllCheckBoxes(self, event):
		for option in self.choiceList.Checked:
			self.choiceList.Check(option, False)


	def selectAllCheckBoxes(self, event):#hardkodirano je promijeni
		toCheck = list(xrange(len(self.checkListBoxInitialItems)))
		self.choiceList.SetChecked(toCheck)


	def generateView(self, event):
		self.testedAssemblers = []
		self.xScale = self.selectionMap[self.xAxisScale.GetSelection()]
		self.yScale = self.selectionMap[self.yAxisScale.GetSelection()]
		self.handleCheckListBox()

		if self.summaryFiles:
			self.loadSummaries()
			try:
				self.generatePlot()
			except Exception:
				dial = wx.MessageDialog(None, "Something went wrong, most probably you are trying to use log scale on negative values, check your input and try to plot again.", 'Plotting error', wx.ICON_ERROR)
				dial.ShowModal()
		else:
			dial = wx.MessageDialog(None, "Load some data before applying changes...", 'Info', wx.ICON_INFORMATION)
			dial.ShowModal()
			
	##############################################################
	#                 Handlers for main menu
	##############################################################
	
	#X Axis basics ----------------------------

	def updateXAxisParameterToRSS(self, event):
		self.xAttribute = 'totalRSS'
		self.updateBasicStatsView()

	def updateXAxisParameterToPSS(self, event):
		self.xAttribute = 'totalPSS'
		self.updateBasicStatsView()

	def updateXAxisParameterToTime(self, event):
		self.xAttribute = 'totalRealTime'
		self.updateBasicStatsView()

	def updateXAxisParameterToVm(self, event):
		self.xAttribute = 'totalVmSize'
		self.updateBasicStatsView()

	def updateXAxisParameterToCPU(self, event):
		self.xAttribute = 'totalCpuTime'
		self.updateBasicStatsView()

	#Y Axis basics ---------------------------

	def updateYAxisParameterToRSS(self, event):
		self.yAttribute = 'totalRSS'
		self.updateBasicStatsView()

	def updateYAxisParameterToPSS(self, event):
		self.yAttribute = 'totalPSS'
		self.updateBasicStatsView()

	def updateYAxisParameterToTime(self, event):
		self.yAttribute = 'totalRealTime'
		self.updateBasicStatsView()

	def updateYAxisParameterToVm(self, event):
		self.yAttribute = 'totalVmSize'
		self.updateBasicStatsView()

	def updateYAxisParameterToCPU(self, event):
		self.yAttribute = 'totalCpuTime'
		self.updateBasicStatsView()
		
	def updateYAxisParameterToCReferenceCoverage(self, event):
		self.yAttribute = 'cReferenceCoverage'
		self.updateBasicStatsView()

	#X Axis wgsim parameters ----------------

	def updateXAxisParameterToL(self, event):
		self.xAttribute = 'l'
		self.updateBasicStatsView()

	def updateXAxisParameterToCov(self, event):
		self.xAttribute = 'cov'
		self.updateBasicStatsView()

	def updateXAxisParameterToN(self, event):
		self.xAttribute = 'N'
		self.updateBasicStatsView()

	def updateXAxisParameterTod(self, event):
		self.xAttribute = 'd'
		self.updateBasicStatsView()

	def updateXAxisParameterToE(self, event):
		self.xAttribute = 'e'
		self.updateBasicStatsView()

	def updateXAxisParameterTor(self, event):
		self.xAttribute = 'r'
		self.updateBasicStatsView()

	def updateXAxisParameterToR(self, event):
		self.xAttribute = 'R'
		self.updateBasicStatsView()

	def updateXAxisParameterToX(self, event):
		self.xAttribute = 'X'
		self.updateBasicStatsView()

	def updateXAxisParameterToA(self, event):
		self.xAttribute = 'A'
		self.updateBasicStatsView()

	def updateXAxisParameterToD(self, event):
		self.xAttribute = 'D'
		self.updateBasicStatsView()
		print 'D'

	#Y - Axis wgsim parameters --------------

	def updateYAxisParameterToL(self, event):
		self.yAttribute = 'l'
		self.updateBasicStatsView()

	def updateYAxisParameterToCov(self, event):
		self.yAttribute = 'cov'
		self.updateBasicStatsView()

	def updateYAxisParameterToN(self, event):
		self.yAttribute = 'N'
		self.updateBasicStatsView()

	def updateYAxisParameterTod(self, event):
		self.yAttribute = 'd'
		self.updateBasicStatsView()

	def updateYAxisParameterToE(self, event):
		self.yAttribute = 'e'
		self.updateBasicStatsView()

	def updateYAxisParameterTor(self, event):
		self.yAttribute = 'r'
		self.updateBasicStatsView()

	def updateYAxisParameterToR(self, event):
		self.yAttribute = 'R'
		self.updateBasicStatsView()

	def updateYAxisParameterToX(self, event):
		self.yAttribute = 'X'
		self.updateBasicStatsView()

	def updateYAxisParameterToA(self, event):
		self.yAttribute = 'A'
		self.updateBasicStatsView()

	def updateYAxisParameterToD(self, event):
		self.yAttribute = 'D'
		self.updateBasicStatsView()

	#X - Axis mixed stats -------------------

	def updateXAxisParameterToAvgs(self, event):
		self.xAttribute = 'mAvgs'
		self.updateBasicStatsView()

	def updateXAxisParameterTo50s(self, event):
		self.xAttribute = 'mN50s'
		self.updateBasicStatsView()

	def updateXAxisParameterToNums(self, event):
		self.xAttribute = 'mNums'
		self.updateBasicStatsView()

	def updateXAxisParameterToLens(self, event):
		self.xAttribute = 'mLens'
		self.updateBasicStatsView()

	def updateXAxisParameterToMaxs(self, event):
		self.xAttribute = 'mMaxs'
		self.updateBasicStatsView()

	#Y - Axis mixed stats ---------------------

	def updateYAxisParameterToAvgs(self, event):
		self.yAttribute = 'mAvgs'
		self.updateBasicStatsView()

	def updateYAxisParameterTo50s(self, event):
		self.yAttribute = 'mN50s'
		self.updateBasicStatsView()

	def updateYAxisParameterToNums(self, event):
		self.yAttribute = 'mNums'
		self.updateBasicStatsView()

	def updateYAxisParameterToLens(self, event):
		self.yAttribute = 'mLens'
		self.updateBasicStatsView()

	def updateYAxisParameterToMaxs(self, event):
		self.yAttribute = 'mMaxs'
		self.updateBasicStatsView()


	#X - Axis scafflod stats ------------------

	def updateXAxisParameterToSTotalNum(self, event):
		self.xAttribute = 'sTotalNum'
		self.updateBasicStatsView()
		
	def updateXAxisParameterToSMaxSize(self, event):
		self.xAttribute = 'sMaxSize'
		self.updateBasicStatsView()

	def updateXAxisParameterToSMinSize(self, event):
		self.xAttribute = 'sMinSize'
		self.updateBasicStatsView()

	def updateXAxisParameterToSAvgSize(self, event):
		self.xAttribute = 'sAvgSize'
		self.updateBasicStatsView()

	def updateXAxisParameterToSMedSize(self, event):
		self.xAttribute = 'sMedSize'
		self.updateBasicStatsView()

	def updateXAxisParameterToSCertainNum(self, event):
		self.xAttribute = 'sCertainNum'
		self.updateBasicStatsView()

	def updateXAxisParameterToSN25(self, event):
		self.xAttribute = 'sN25'
		self.updateBasicStatsView()

	def updateXAxisParameterToSN50(self, event):
		self.xAttribute = 'sN50'
		self.updateBasicStatsView()

	def updateXAxisParameterToSN56(self, event):
		self.xAttribute = 'sN56'
		self.updateBasicStatsView()

	def updateXAxisParameterToSN75(self, event):
		self.xAttribute = 'sN75'
		self.updateBasicStatsView()

	def updateXAxisParameterToSEM(self, event):
		self.xAttribute = 'sEM'
		self.updateBasicStatsView()

	def updateXAxisParameterToSEN(self, event):
		self.xAttribute = 'sEN'
		self.updateBasicStatsView()

	def updateXAxisParameterToSQual(self, event):
		self.xAttribute = 'sQual'
		self.updateBasicStatsView()

	#Y - Axis scaffold stats

	
	def updateYAxisParameterToSTotalNum(self, event):
		self.yAttribute = 'sTotalNum'
		self.updateBasicStatsView()
		
	def updateYAxisParameterToSMaxSize(self, event):
		self.yAttribute = 'sMaxSize'
		self.updateBasicStatsView()

	def updateYAxisParameterToSMinSize(self, event):
		self.yAttribute = 'sMinSize'
		self.updateBasicStatsView()

	def updateYAxisParameterToSAvgSize(self, event):
		self.yAttribute = 'sAvgSize'
		self.updateBasicStatsView()

	def updateYAxisParameterToSMedSize(self, event):
		self.yAttribute = 'sMedSize'
		self.updateBasicStatsView()

	def updateYAxisParameterToSCertainNum(self, event):
		self.yAttribute = 'sCertainNum'
		self.updateBasicStatsView()

	def updateYAxisParameterToSN25(self, event):
		self.yAttribute = 'sN25'
		self.updateBasicStatsView()

	def updateYAxisParameterToSN50(self, event):
		self.yAttribute = 'sN50'
		self.updateBasicStatsView()

	def updateYAxisParameterToSN56(self, event):
		self.yAttribute = 'sN56'
		self.updateBasicStatsView()

	def updateYAxisParameterToSN75(self, event):
		self.yAttribute = 'sN75'
		self.updateBasicStatsView()

	def updateYAxisParameterToSEM(self, event):
		self.yAttribute = 'sEM'
		self.updateBasicStatsView()

	def updateYAxisParameterToSEN(self, event):
		self.yAttribute = 'sEN'
		self.updateBasicStatsView()

	def updateYAxisParameterToSQual(self, event):
		self.yAttribute = 'sQual'
		self.updateBasicStatsView()

	#X - Axis contig stats ---------------------

	def updateXAxisParameterToCTotalNum(self, event):
		self.xAttribute = 'cTotalNum'
		self.updateBasicStatsView()

	def updateXAxisParameterToCBiggerThen(self, event):
		self.xAttribute = 'cBiggerThen'
		self.updateBasicStatsView()

	def updateXAxisParameterToCTotalLen(self, event):
		self.xAttribute = 'cTotalLen'
		self.updateBasicStatsView()

	def updateXAxisParameterToCMaxLen(self, event):
		self.xAttribute = 'cMaxLen'
		self.updateBasicStatsView()
		
	def updateXAxisParameterToCMinLen(self, event):
		self.xAttribute = 'cMinLen'
		self.updateBasicStatsView()

	def updateXAxisParameterToCAvgLen(self, event):
		self.xAttribute = 'cAvgLen'
		self.updateBasicStatsView()

	def updateXAxisParameterToCMedLen(self, event):
		self.xAttribute = 'cMedLen'
		self.updateBasicStatsView()

	def updateXAxisParameterToCESize(self, event):
		self.xAttribute = 'cESize'
		self.updateBasicStatsView()

	def updateXAxisParameterToCN25(self, event):
		self.xAttribute = 'cN25'
		self.updateBasicStatsView()

	def updateXAxisParameterToCN50(self, event):
		self.xAttribute = 'cN50'
		self.updateBasicStatsView()

	def updateXAxisParameterToCN56(self, event):
		self.xAttribute = 'cN56'
		self.updateBasicStatsView()

	def updateXAxisParameterToCN75(self, event):
		self.xAttribute = 'cN75'
		self.updateBasicStatsView()
		print 'cN75'

	#Y - Axis contig stats ----------------------

	def updateYAxisParameterToCTotalNum(self, event):
		self.yAttribute = 'cTotalNum'
		self.updateBasicStatsView()

	def updateYAxisParameterToCBiggerThen(self, event):
		self.yAttribute = 'cBiggerThen'
		self.updateBasicStatsView()

	def updateYAxisParameterToCTotalLen(self, event):
		self.yAttribute = 'cTotalLen'
		self.updateBasicStatsView()

	def updateYAxisParameterToCMaxLen(self, event):
		self.yAttribute = 'cMaxLen'
		self.updateBasicStatsView()
		
	def updateYAxisParameterToCMinLen(self, event):
		self.yAttribute = 'cMinLen'
		self.updateBasicStatsView()

	def updateYAxisParameterToCAvgLen(self, event):
		self.yAttribute = 'cAvgLen'
		self.updateBasicStatsView()

	def updateYAxisParameterToCMedLen(self, event):
		self.yAttribute = 'cMedLen'
		self.updateBasicStatsView()

	def updateYAxisParameterToCESize(self, event):
		self.yAttribute = 'cESize'
		self.updateBasicStatsView()

	def updateYAxisParameterToCN25(self, event):
		self.yAttribute = 'cN25'
		self.updateBasicStatsView()

	def updateYAxisParameterToCN50(self, event):
		self.yAttribute = 'cN50'
		self.updateBasicStatsView()

	def updateYAxisParameterToCN56(self, event):
		self.yAttribute = 'cN56'
		self.updateBasicStatsView()

	def updateYAxisParameterToCN75(self, event):
		self.yAttribute = 'cN75'
		self.updateBasicStatsView()

	def loadSummaries(self): #from this summary parsers generate plot
		self.summaryParsers = []
		for summaryFile in self.activeSummaryFiles:
			summaryData = summaryParser.SummaryParser(summaryFile)
			self.summaryParsers.append(summaryData)

	##############################################################
	#                     Export functions
	##############################################################

	def exportBasicStatsToPDF(self, event):
		print "[AN:] Exporting basic statistics to pdf..."
		if self.plotDirectory == '':
			dial = wx.MessageDialog(None, "Directory for storing plots is not specified!", 'Error', wx.ICON_ERROR)
			dial.ShowModal()
			print "[AN:] Aborted."
		else:
			fileName = self.xAttribute + "_" + self.yAttribute + "_bPlot.pdf"
			path = os.path.join(self.plotDir, fileName)
			self.figure.savefig(path, dpi=300)
			dial = wx.MessageDialog(None, "Plot exported successfully...", 'Info', wx.ICON_INFORMATION)
			dial.ShowModal()
			print "[AN:] Done..."

	def exportBasicStatsToPng(self, event):
		print "[AN:] Exporting basic statistics to png..."
		if self.plotDir == '':
			dial = wx.MessageDialog(None, "Directory for storing plots is not specified!", 'Error', wx.ICON_ERROR)
			dial.ShowModal()
			print "[AN:] Aborted."
		else:
			fileName = self.xAttribute + "_" + self.yAttribute + "_bPlot.png"
			path = os.path.join(self.plotDir, fileName)
			self.figure.savefig(path, dpi=300)
			dial = wx.MessageDialog(None, "Plot exported successfully...", 'Info', wx.ICON_INFORMATION)
			dial.ShowModal()
			print "[AN:] Done..."

	def exportMappedReadsToPDF(self, event):
		print "[AN:] Exporting mapped reads to pdf..."
		if self.readFile == '':
			dial = wx.MessageDialog(None, "Input is not specified, cannot export sample plot!", 'Error', wx.ICON_ERROR)
			dial.ShowModal()
			print "[AN:] Aborted."
		else:
			fileName = self.readFile.split('\\')[-1] + "_" + self.readPlotName.GetLabel() + "_mapping.pdf" #pazi linux
			path = os.path.join(self.plotDir, fileName)
			self.readFigure.savefig(path, dpi=300)
			dial = wx.MessageDialog(None, "Plot exported successfully...", 'Info', wx.ICON_INFORMATION)
			dial.ShowModal()
			print "[AN:] Done..."

	def exportMappedReadsToPNG(self, event):
		print "[AN:] Exporting mapped reads to png..."
		if self.readFile == '':
			dial = wx.MessageDialog(None, "Input file is not specified, cannot export sample plot!", 'Error', wx.ICON_ERROR)
			dial.ShowModal()
			print "[AN:] Aborted."
		else:
			fileName = self.readFile.split('\\')[-1] + "_" + self.readPlotName.GetLabel() + "_mapping.png" #pazi linux
			path = os.path.join(self.plotDir, fileName)
			self.readFigure.savefig(path, dpi=300)
			dial = wx.MessageDialog(None, "Plot exported successfully...", 'Info', wx.ICON_INFORMATION)
			dial.ShowModal()
			print "[AN:] Done..."
		

	def setPlotDirectory(self, event):
		self.plotDir = self.plotDirectory.GetPath()
		print self.plotDir


	def setYAxisScale(self, event):
		pass

	def showSomeInfo(self, event):
		message1 = "Click the browse button and select your summary folder.\nCheckboxes will be filled automatically."
		message2 = " Once when the data is loaded the sample plot on the second tab will be replaced with real one."
		dial = wx.MessageDialog(None, message1 + message2, 'Info', wx.OK)
		dial.ShowModal()

	def getMappingHelp(self, event):
		message = "Pick a reference genome file and contig file to start the simulation..."
		dial = wx.MessageDialog(None, message, 'Info', wx.OK)
		dial.ShowModal()
		
	def threadTarget(self, referenceGenomeFile, contigFile):
		self.cmd1 = "nucmer -maxmatch -c 100 -p nucmer "+referenceGenomeFile+" "+contigFile
		self.cmd2 = "show-coords -r -c -l nucmer.delta > nucmer.coords"
		self.cmd3 = "show-snps -C nucmer.delta > nucmer.snps"
		self.cmd4 = "show-tiling nucmer.delta > nucmer.tiling"
		self.cmd5 = "delta-filter -m nucmer.delta > nucmer.delta.m"
		self.cmd6 = "mummerplot nucmer.delta.m"
		self.cmd7 = "dnadiff -d nucmer.delta"
		os.system(self.cmd1)
		os.system(self.cmd2)
		os.system(self.cmd3)
		os.system(self.cmd4)
		os.system(self.cmd5)
		

	def runMummer(self,event):
		referenceGenomeFile = self.referencePicker.GetPath()
		contigFile = self.contigPicker.GetPath()
		if referenceGenomeFile == '' or contigFile == '':
			message = "You haven't selected necessary files..."
			dial = wx.MessageDialog(None, message, 'Info', wx.OK)
			dial.ShowModal()
		else:
		   execute = "python"
		   print referenceGenomeFile, contigFile
		   arg1 = "MUMMERCaller.py"
		   arg2 = referenceGenomeFile
		   arg3 = contigFile
		   subprocess.call([execute, arg1, arg2, arg3])

	


app = wx.App(False)
frame = BenchmarkMainFrame(None)
frame.Show(True)
app.MainLoop()
