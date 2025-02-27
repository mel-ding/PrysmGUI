from os.path import basename
import psm.cellulose.sensor as sensor
import psm.agemodels.banded as banded
import psm.aux_functions.analytical_err_simple as analytical_err_simple
import psm.aux_functions.analytical_error as analytical_error

import tkinter as tk 	
from tkinter.filedialog import asksaveasfilename, askopenfilename
import numpy as np
from scipy.stats.mstats import mquantiles
import nitime.algorithms as tsa
from pandas import DataFrame

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

LARGE_FONT = ("Verdana", 25) 

class Grapher(tk.Frame):

	def __init__(self, parent, *args, **kwargs):
		tk.Frame.__init__(self, parent, *args, **kwargs)
		self.parent = parent
		self.parent.title("PRYSM: Cellulose PSM")

		# Initialize empty arrays for data saving 
		self.time = np.array([])
		self.temp=np.array([])
		self.precip=np.array([])
		self.rh=np.array([])
		self.d180v=np.array([])
		self.d180s=np.array([])
		self.d180p=np.array([])
		self.Xp = np.array([])
		self.ageq1 = np.array([])
		self.gaussq1 = np.array([])
		self.simpleq1 = np.array([])
		self.simpleq2 = np.array([])
		self.prepSuccess = False 	# keep track of successful data prep
		self.PS = False				# keep track of what graph is displaying 
		self.newSensorData = True 	# do we need to calculate new sensor data?
		self.newPSData = True 		# do we need to calculate new PS data?
		self.R1 = -10.0				
		self.R2 = -10.0 				
		self.R3 = -10.0 	
		self.model = -1		
		self.d180 = False 			

		# =========================================================================================
		# SENSOR DATA 
		# =========================================================================================

		# Title of Page 
		label = tk.Label(root, text="PRYSM: Cellulose PSM", font=LARGE_FONT)
		label.grid(sticky="E")

		rowIdx = 1

		# =========================================================================================
		# FLAG TOGGLES 
		# =========================================================================================

		# Rodens or Evans 
		self.v = tk.IntVar()
		tk.Label(root, text = "Pick a model:").grid(row=rowIdx, column=0, sticky="E")
		tk.Radiobutton(root, text="Roden", variable=self.v, value=0).grid(row=rowIdx, column=1, sticky="w")
		rowIdx += 1
		tk.Radiobutton(root, text="Evans", variable=self.v, value=1).grid(row=rowIdx, column=1, sticky="w")
		rowIdx += 1

		# d180 inputs?
		self.inputs = tk.IntVar()
		tk.Label(root, text = "Do you have d180 data?").grid(row=rowIdx, column=0, sticky="E")
		tk.Radiobutton(root, text="No", variable=self.inputs, value=0).grid(row=rowIdx, column=1, sticky="w")
		rowIdx += 1
		tk.Radiobutton(root, text="Yes", variable=self.inputs, value=1).grid(row=rowIdx, column=1, sticky="w")
		rowIdx += 1

		# =========================================================================================
		# ERROR OPTIONS
		# =========================================================================================

		# Age Uncertainties
		self.ageErrorVal = tk.BooleanVar() 
		self.ageErrorVar = tk.Checkbutton(root, variable=self.ageErrorVal)
		self.ageErrorVar.grid(row=rowIdx, column=2)
		tk.Label(root, text = "Add age uncertainties with error rate:").grid(row=rowIdx, column=0, sticky="E")
		self.ageErrorEntry = tk.Entry(root)
		self.ageErrorEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Simple Analytical Errors
		self.simAltErrorVal = tk.BooleanVar() 
		self.simAltErrorVar = tk.Checkbutton(root, variable=self.simAltErrorVal)
		self.simAltErrorVar.grid(row=rowIdx, column=2)
		tk.Label(root, text = "Add simple analytical error with precision:").grid(row=rowIdx, column=0, sticky="E")
		self.simAltErrorEntry = tk.Entry(root)
		self.simAltErrorEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Gaussian Analytical Errors
		self.altErrorVal = tk.BooleanVar() 
		self.altErrorVar = tk.Checkbutton(root, variable=self.altErrorVal)
		self.altErrorVar.grid(row=rowIdx, column=2)
		tk.Label(root, text = "Add gaussian analytical error with precision:").grid(row=rowIdx, column=0, sticky="E")
		self.altErrorEntry = tk.Entry(root)
		self.altErrorEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# =========================================================================================
		# BUTTONS 
		# =========================================================================================

		# Allows user to upload data. 
		tk.Label(root, 
			text="Upload a single CSV file with columns and the headers \n \"TIME\", \"PRECIPITATION\", \"HUMIDITY\", \"TEMPERATURE\", \n \"D180_VAPOR\", \"D180_SOIL\", \"D180_PRECIPITATION\"."
			).grid(row=rowIdx, columnspan=3, rowspan=4)
		rowIdx += 4
		tk.Label(root, text = "Click to upload your data:").grid(row=rowIdx, column=0, sticky="E")
		graphButton = tk.Button(root, text="Upload Data", command=self.uploadData)
		graphButton.grid(row=rowIdx, column=1, ipadx=30, ipady=3, sticky="W")
		rowIdx += 1

		# Shows the name of the current uploaded file, if any. 
		tk.Label(root, text="Current File Uploaded:").grid(row=rowIdx, column=0, sticky="E")
		self.currentFileLabel = tk.Label(root, text="No file")
		self.currentFileLabel.grid(row=rowIdx, column=1, columnspan=2, sticky="W")
		rowIdx+=1

		# Creates graph with given inputs
		tk.Label(root, text = "Click to create your own graph:").grid(row=rowIdx, column=0, sticky="E")
		graphButton = tk.Button(root, text="Generate Graph", command=self.generateGraph)
		graphButton.grid(row=rowIdx, column=1, ipadx=20, ipady=3, sticky="W")
		rowIdx+=1

		# Clears graph
		tk.Label(root, text = "Click to clear the graph:").grid(row=rowIdx, column=0, sticky="E")
		clearGraphButton = tk.Button(root, text="Clear Graph", command=self.clearGraph)
		clearGraphButton.grid(row=rowIdx, column=1, ipadx=33, ipady=3, sticky="W")
		rowIdx+=1

		# Clears entries
		tk.Label(root, text = "Click to clear the entries:").grid(row=rowIdx, column=0, sticky="E")
		clearEntryButton = tk.Button(root, text="Clear Entries", command=self.clearEntries)
		clearEntryButton.grid(row=rowIdx, column=1, ipadx=30, ipady=3, sticky="W")
		rowIdx+=1 

		# =========================================================================================
		# GRAPH
		# =========================================================================================
		PSGraphButton = tk.Button(root, text="Power Spectrum", command=self.generatePS)
		PSGraphButton.grid(row=0, column=9, ipadx=20, ipady=3, sticky="S")

		SensorGraphButton = tk.Button(root, text="Sensor", command=self.generateGraph)
		SensorGraphButton.grid(row=0, column=11, ipadx=20, ipady=3, sticky="S")

		self.f = Figure(figsize=(10,5), dpi=100)
		self.plt = self.f.add_subplot(111)
		self.canvas = FigureCanvasTkAgg(self.f, root)
		self.canvas.get_tk_widget().grid(row=1, column=3, rowspan=16, columnspan=15, sticky="nw")
		self.plt.set_title(r'SENSOR',fontsize=12)
		self.plt.set_xlabel('Time')
		self.plt.set_ylabel('Simulated Cellulose Data')

		# =========================================================================================
		# SAVE OPTIONS 
		# =========================================================================================	
		saveRowIdx = 19
		tk.Label(root, text = "Save sensor data as:").grid(row=saveRowIdx, column=5, sticky="E")
		dataTxtbutton = tk.Button(root, text=".txt", command=self.saveTxtData)
		dataTxtbutton.grid(row=saveRowIdx, column=7, ipadx=20, ipady=3, sticky="W")
		dataCsvbutton = tk.Button(root, text=".csv", command=self.saveCsvData)
		dataCsvbutton.grid(row=saveRowIdx, column=8, ipadx=20, ipady=3, sticky="W")

		tk.Label(root, text = "Save uncertainty data as:").grid(row=saveRowIdx, column=9, sticky="E")
		errorTxtbutton = tk.Button(root, text=".txt", command=self.saveTxtErrors)
		errorTxtbutton.grid(row=saveRowIdx, column=11, ipadx=20, ipady=3, sticky="W")
		errorCsvbutton = tk.Button(root, text=".csv", command=self.saveCsvErrors)
		errorCsvbutton.grid(row=saveRowIdx, column=12, ipadx=20, ipady=3, sticky="W")
		saveRowIdx += 1

		tk.Label(root, text = "Save power spectrum data as:").grid(row=saveRowIdx, column=5, sticky="E")
		errorTxtbutton = tk.Button(root, text=".txt", command=self.saveTxtPSData)
		errorTxtbutton.grid(row=saveRowIdx, column=7, ipadx=20, ipady=3, sticky="W")
		errorCsvbutton = tk.Button(root, text=".csv", command=self.saveCsvPSData)
		errorCsvbutton.grid(row=saveRowIdx, column=8, ipadx=20, ipady=3, sticky="W")

		tk.Label(root, text = "Save graph as:").grid(row=saveRowIdx, column=9, sticky="E")
		graphPNGbutton = tk.Button(root, text=".png", command=self.savePngGraph)
		graphPNGbutton.grid(row=saveRowIdx, column=11, ipadx=20, ipady=3, sticky="W")
		graphPDFbutton = tk.Button(root, text=".pdf", command=self.savePdfGraph)
		graphPDFbutton.grid(row=saveRowIdx, column=12, ipadx=20, ipady=3, sticky="W")
		saveRowIdx += 1

	"""
	Saves power spectrum data as text file. 
	"""
	def saveTxtPSData(self): 
		if (self.Cd18Of.size != 0):
			file1 = asksaveasfilename(initialfile="Cd180f", defaultextension=".txt")
			if file1:
				np.savetxt(file1, self.Cd18Of, newline=" ")
				tk.messagebox.showinfo("Sucess", "Saved power spectrum frequency data")
		if (self.Cd18Opsd_mt.size != 0):
			file2 = asksaveasfilename(initialfile="Cd18Opsd_mt", defaultextension=".txt")
			if file2:
				np.savetxt(file2, self.Cd18Opsd_mt, newline=" ")
				tk.messagebox.showinfo("Sucess", "Saved power spectral density data")

	"""
	Saves power spectrum data as csv file. 
	"""
	def saveCsvPSData(self): 
		df = DataFrame({})
		if (self.Cd18Of.size != 0):
			df["Cd180f"] = self.Cd18Of
		if (self.Cd18Opsd_mt.size != 0):
			df["Cd18Opsd_mt"] = self.Cd18Opsd_mt
		file = asksaveasfilename(initialfile="PowerSpectrumData.csv", defaultextension=".csv")
		if file:
			df.to_csv(file, index=False)
			tk.messagebox.showinfo("Sucess", "Saved power spectrum data")

	"""
	Saves current graph as pdf.
	"""
	def savePdfGraph(self):
		file = asksaveasfilename(initialfile="Figure.pdf", defaultextension=".pdf")
		if file:
			self.f.savefig(file)
			tk.messagebox.showinfo("Sucess", "Saved graph")

	"""
	Saves current graph as png.
	"""
	def savePngGraph(self):
		file = asksaveasfilename(initialfile="Figure.png", defaultextension=".png")
		if file:
			self.f.savefig(file)
			tk.messagebox.showinfo("Sucess", "Saved graph")

	"""
	Saves sensor data into a text file. 
	"""
	def saveTxtData(self):
		file = asksaveasfilename(initialfile="SensorData.txt", defaultextension=".txt")
		if file:
			np.savetxt(file, self.cell, newline=" ")
			tk.messagebox.showinfo("Sucess", "Saved simulated data")
	
	"""
	Saves sensor data into a csv file. 
	"""
	def saveCsvData(self):
		df = DataFrame({"Sensor": self.cell})
		file = asksaveasfilename(initialfile="SensorData.csv", defaultextension=".csv")
		if file:
			df.to_csv(file, index=False)
			tk.messagebox.showinfo("Sucess", "Saved uncertainty data")

	"""
	Saves error data into a text file. 
	"""
	def saveCsvErrors(self):
		df = DataFrame({})
		if (self.ageq1.size != 0):
			df["Age_Q1"] = self.ageq1[:,0]
			df["Age_Q2"] = self.ageq1[:,1]

		if (self.gaussq1.size != 0):
			df["GaussianAnalytical_Q1"] = self.gaussq1[:,0]
			df["GaussianAnalytical_Q2"] = self.gaussq1[:,1]
		
		if (self.simpleq1.size != 0 and self.simpleq2.size != 0):
			df["SimpleAnalytical_Q1"] = self.simpleq1
			df["SimpleAnalytical_Q2"] = self.simpleq2
			
		file = asksaveasfilename(initialfile="uncertainties.csv", defaultextension=".csv")
		if file:
			df.to_csv(file, index=False)
			tk.messagebox.showinfo("Sucess", "Saved uncertainty data")

	"""
	Saves error data into a text file. 
	"""
	def saveTxtErrors(self):	
		if (self.ageq1.size != 0):
			file1 = asksaveasfilename(initialfile="age_uncertainties", defaultextension=".txt")
			if file1:
				np.savetxt(file1, (self.ageq1[:,0], self.ageq1[:,1]), newline=" ")
				tk.messagebox.showinfo("Sucess", 
					"Saved age error data where Q1 is the first row and Q2 is the second row.")
		if (self.gaussq1.size != 0):
			file2 = asksaveasfilename(initialfile="gaussian_analytical_uncertainties", defaultextension=".txt")
			if file2:			
				np.savetxt(file2, (self.gaussq1[:,0], self.gaussq1[:,1]), newline=" ")
				tk.messagebox.showinfo("Sucess", 
					"Saved gaussian analytical error data where Q1 is the first row and Q2 is the second row.")
		if (self.simpleq1.size != 0 and self.simpleq2.size != 0):
			file3 = asksaveasfilename(initialfile="simple_analytical_uncertainties", defaultextension=".txt")
			if file3:	
				np.savetxt(file3, (self.simpleq1, self.simpleq2), newline=" ")
				tk.messagebox.showinfo("Sucess", 
					"Saved simple analytical error data where Q1 is the first row and Q2 is the second row.")

	"""
	Takes a CSV file and saves the Time, SSS, and SST values. 
	"""
	def uploadData(self):
		# Open the file choosen by the user 
		filename = askopenfilename(filetypes = (("csv files","*.csv"),))
		data = np.genfromtxt(filename, delimiter = ",", names=True, dtype=None)
		self.currentFileLabel.configure(text=basename(filename))
		# Get the entry fields.  		
		self.time=data['TIME']
		self.temp=data['TEMPERATURE']
		self.precip=data['PRECIPITATION']
		self.rh=data['HUMIDITY']
		self.d180v=data['D180_VAPOR']
		self.d180s=data['D180_SOIL']
		self.d180p=data['D180_PRECIPITATION']
		self.newData = True	
		if (self.time.shape != self.precip.shape != self.rh.shape != self.d180s.shape != self.d180v.shape != self.d180p.shape):
			tk.messagebox.showerror("Error", "Invalid Data: Data inputs are different lengths.")	
			return
	
	"""
	Take the input data and calculate the sensor data. 
	"""	
	def dataPrep(self): 

		# Make sure user has uploaded data 
		if (self.time.size == 0 or self.precip.size == 0 or self.rh.size == 0 or self.temp.size == 0):
			tk.messagebox.showerror("Error", "Missing input data")
			self.prepSuccess = False
			return

		# If user uploaded d180 data, make sure valid
		if (self.inputs.get() == 1):
			if (self.d180s.size == 0 or self.d180v.size == 0 or self.d180p.size == 0): 
				tk.messagebox.showerror("Error", "Missing input data")
				self.prepSuccess = False
				return

		# Clear whatever is currently on the canvas 
		self.plt.clear()

		# Get which model user selected
		modelInput = self.v.get()

		# Get indication is user input d180 data or not 
		d180_data = False
		if (self.inputs.get() == 1):
			d180_data = True 

		if (self.model != modelInput or self.d180 != d180_data): 
			self.newSensorData = False
			self.newPSData = True
			self.model = modelInput 
			self.d180 = d180_data 
			# Fill cellulose array with data same size as input vectors.
			self.cell = sensor.cellulose_sensor(self.time, self.temp, self.precip, self.rh, 
				self.d180s, self.d180p, self.d180v, flag=self.model, iso=self.d180)
			# Reshape cellulose data for uncertainty calculations 
			self.X = self.cell
			self.X = self.X.reshape(len(self.X),1)
		
		self.prepSuccess = True

	"""
	Generates coral data based on input data and model, and graphs result.
	"""
	def generateGraph(self): 

		self.dataPrep()

		# if data prep was unsucessful, do not try to generate graph.
		if (self.prepSuccess == False): 
			return

		# change from 3 subplots to 1 
		if (self.PS == True): 
			self.ax.remove()
			self.ax2.remove()
			self.PS = False 

		# Plot age uncertainties - if selected 
		if (self.ageErrorVal.get()):
			# Get input error rate
			rateRaw = self.ageErrorEntry.get()
			if rateRaw == "":
				tk.messagebox.showerror("Error", "Invalid Rate Value")
				return
			# if first time selecting rate, remember 
			if (self.R1 == -10):
				rate = float(rateRaw)
				self.R1 = rate 				
			# if rate changed or new data, recalculate 
			if (self.R1 != float(rateRaw) or self.newSensorData == True):
				rate = float(rateRaw)
				self.R1 = rate 
				# Calculate the age uncertanties
				tp, self.Xp, tmc = banded.bam_simul_perturb(self.X, self.time, param=[rate, rate])
				self.ageq1=mquantiles(self.Xp, prob=[0.025,0.975],axis=1)
				q2=self.time
				# Graph quantiles 
				self.plt.fill_between(q2, self.ageq1[:,0], self.ageq1[:,1],
					label='1000 Age-Perturbed Realizations, CI', facecolor='gray',alpha=0.5)

		# Plot simple analytical errors - if selected
		if (self.simAltErrorVal.get()):
			# Get input error rate
			sigmaRaw = self.simAltErrorEntry.get()
			if sigmaRaw == "":
				tk.messagebox.showerror("Error", "Invalid Rate Value")
				return
			# if first time selecting rate, remember 
			if (self.R2 == -10):
				sigma = float(sigmaRaw)
				self.R2 = sigma 				
			# if rate changed or new data, recalculate 
			if (self.R2 != float(sigmaRaw) or self.newSensorData == True):
				sigma = float(sigmaRaw)
				self.R2 = sigma 
				self.simpleq1, self.simpleq2 = analytical_err_simple.analytical_err_simple(self.X,sigma)
				# Reshape quanties for graphing
				self.simpleq1 = self.simpleq1.reshape(len(self.time))
				self.simpleq2 = self.simpleq2.reshape(len(self.time))
				# Graph quantiles 
				self.plt.fill_between(self.time, self.simpleq1, self.simpleq2, 
					label='100 Analytical Error Realizations, CI', facecolor='darkgray',alpha=0.5)

		# Plot gaussian analytical errors - if selected
		if (self.altErrorVal.get()):
			# Get input error rate
			sigmaRaw = self.altErrorEntry.get()
			if sigmaRaw == "":
				tk.messagebox.showerror("Error", "Invalid Rate Value")
				return
			# if first time selecting rate, remember 
			if (self.R3 == -10):
				sigma = float(sigmaRaw)
				self.R3 = sigma 				
			# if rate changed or new data, recalculate 
			if (self.R3 != float(sigmaRaw) or self.newSensorData == True):
				sigma = float(sigmaRaw)
				self.R3 = sigma 
				inputs = len(self.X)
				self.X = self.cell
				self.X = self.X.reshape(inputs,1)
				Xn = analytical_error.analytical_error(self.X, sigma, inputs)
				Xn = Xn[:,0,:].reshape(inputs, inputs) 
				self.gaussq1=mquantiles(Xn,prob=[0.025,0.975],axis=1)
				q2=self.time
				self.plt.fill_between(q2,self.gaussq1[:,0],self.gaussq1[:,1],
					label='100 Gaussian Analytical Error Realizations, CI', facecolor='darkgray',alpha=0.5)

		# Plot the graph	
		self.plt = self.f.add_subplot(111)
		self.plt.set_title(r'SENSOR')
		self.plt.set_xlabel('Time')
		self.plt.set_ylabel('Simulated Cellulose Data')
		self.plt.plot(self.time, self.cell, color='DarkGreen')
		self.canvas = FigureCanvasTkAgg(self.f, root)
		self.canvas.get_tk_widget().grid(row=1, column=3, rowspan=16, columnspan=15, sticky="nw")	
		self.canvas.draw()

	"""
	Generates power spectrum of coral model. 
	"""
	def generatePS(self):

		self.dataPrep()

		# if data prep was unsucessful, do not try to generate graph.
		if (self.prepSuccess == False): 
			return

		# clear whatever is currently on the canvas 
		self.plt.clear()
		self.f.clear()
		self.PS = True 

		if (self.newPSData == True):	
			self.TC = self.cell - np.mean(self.cell)
			self.TCf, self.TCpsd_mt, self.TCnu = tsa.multi_taper_psd(self.TC, Fs=1.0, 
				adaptive=False, jackknife=False)

		self.t = np.arange(1000,2005,1)
		self.dt = 1.0
		
		# Get input error rate, if specified 
		rateRaw = self.ageErrorEntry.get()
		if rateRaw == "":
			rate = 0.02
		else:
			rate = float(rateRaw)

		# if rate hasn't changed and already calculated, no need to recalculate 
		if (self.R1 != rate or self.newPSData == True):
			self.R1 = rate
			self.tp, self.Xp, self.tmc = banded.bam_simul_perturb(self.X, self.time, param=[rate,rate], 
				name='poisson', ns=100, resize=0)
			self.Xpm = self.Xp - np.mean(self.Xp, axis=0)
			# DO SAME CALCULATION IN LOOP FOR ALL AGE UNCERTAINTY VECTORS
			self.Xpf=np.zeros((len(self.TCf), len(self.Xpm[1])))
			self.Xpsd_mt=np.zeros((len(self.TCf), len(self.Xpm[1])))
			self.Xpnu=np.zeros((len(self.TCf), len(self.Xpm[1])))
			for i in range(len(self.Xpm[1])):
				self.Xpf[:,i], self.Xpsd_mt[:,i], self.Xpnu[:,i] = tsa.multi_taper_psd(self.Xpm[:,i], Fs=1.0,
					adaptive=False, jackknife=False)
			# Compute quantiles for spectra
			self.q1 = mquantiles(self.Xpsd_mt, prob=[0.025,0.975], axis=1)
			# x axis for quantile of spectra:
			self.q2 = self.Xpf[:,0]		
			self.newPSData = False

		# =========================================================================================
		
		self.ax=self.f.add_subplot(211)
		self.ax.spines["top"].set_visible(False)  
		self.ax.spines["right"].set_visible(False) 
		self.ax.loglog(self.TCf, self.TCpsd_mt, color='DarkGreen') 

		# Change the x-tick labels to something reasonable
		pertick = [500,200,100,50,20,10,8,6,4,2]
		xtick = [500,200,100,50,20,10,8,6,4,2]
		for i in range(len(pertick)):
			xtick[i] = 1.0/pertick[i]
		pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])
		
		self.ax.set_title(r'ENVIRONMENT', color='gray')
		self.ax.set_xlabel(r'Period (Years)')
		self.ax.set_ylabel(r'PSD')
		self.ax.minorticks_off()
		self.ax.set_xticks(xtick)
		self.ax.set_xticklabels(pertick_labels)
		self.ax.grid('on',axis='y',color='DimGray')
		self.ax.set_xlim([1./550.,0.4])
		self.ax.set_ylim([1e-2, 1e0])
		self.ax.legend(loc=3,fontsize=11,frameon=False)

		# =========================================================================================

		self.ax2=self.f.add_subplot(212)
		self.ax2.spines["top"].set_visible(False)  
		self.ax2.spines["right"].set_visible(False)  
		# Observation purturbed age ensemble here.
		self.ax2.fill_between(self.q2, self.q1[:,0], self.q1[:,1], label='1000 Age-Perturbed Realizations, CI',
			facecolor='gray', alpha=0.5)
		self.ax2.loglog(self.TCf, self.TCpsd_mt, color='DarkGreen')

		self.ax2.set_title(r'OBSERVATION', fontsize=11, color='gray')
		self.ax2.set_xlabel(r'Frequency (Years)')
		self.ax2.set_ylabel(r'PSD')
		self.ax2.minorticks_off()
		self.ax2.grid('on',axis='y', color='DimGray')
		self.ax.set_xticks(xtick)
		self.ax.set_xticklabels(pertick_labels)  
		self.ax2.set_xlim([1./550., 0.4])
		self.ax2.set_ylim([1e-2, 1e0])
		self.ax2.legend(loc=3, fontsize=11, frameon=False)

		# =========================================================================================

		self.f.subplots_adjust(hspace=.95)
		self.canvas = FigureCanvasTkAgg(self.f, root)
		self.canvas.get_tk_widget().grid(row=1, column=3, rowspan=20, columnspan=15, sticky="nw")

	"""
	Clears all content in the graph. 
	"""
	def clearGraph(self):
		# clear the axes from power spectrum graph 
		if (self.PS == True):
			self.ax.cla()
			self.ax2.cla() 
		self.plt.clear()
		self.canvas.draw()

	"""
	Clears entries the user has input. 
	"""
	def clearEntries(self):
		self.ageErrorVal.set(0)
		self.simAltErrorVal.set(0)
		self.altErrorVal.set(0)
		self.simAltErrorEntry.delete(0, tk.END)
		self.ageErrorEntry.delete(0, tk.END)
		self.altErrorEntry.delete(0, tk.END)
		self.v.set(None)
		self.inputs.set(None)


"""
Initialize Application.
"""
if __name__ == "__main__":
	root = tk.Tk()
	Grapher(root)
	root.geometry("1600x700+0+50")
	root.mainloop()

