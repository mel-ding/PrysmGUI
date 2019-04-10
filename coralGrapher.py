import psm.coral.sensor as sensor
import psm.agemodels.banded as banded
import psm.aux_functions.analytical_err_simple as analytical_err_simple
import psm.aux_functions.analytical_error as analytical_error

import tkinter as tk 	
from tkinter import ttk, filedialog
import numpy as np
from scipy.stats.mstats import mquantiles
import nitime.algorithms as tsa
import pandas as pd

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

LARGE_FONT = ("Verdana", 25) 

class Grapher(tk.Frame):

	def __init__(self, parent, *args, **kwargs):
		tk.Frame.__init__(self, parent, *args, **kwargs)
		self.parent = parent
		self.parent.title("Coral PRSYM Model")

		# Initialize empty arrays for data saving 
		self.time = np.array([])
		self.SST = np.array([])
		self.SSS = np.array([])
		self.Xp = np.array([])
		self.ageq1 = np.array([])
		self.gaussq1 = np.array([])
		self.simpleq1 = np.array([])
		self.simpleq2 = np.array([])
		self.example = False 		# keep track of example data 
		self.PS = False				# keep track of what graph is displaying 
		self.newData = False 		# did the user upload new data?
		self.R1 = -10.0				
		self.R2 = -10.0 				
		self.R3 = -10.0 				

		# =========================================================================================
		# SENSOR DATA 
		# =========================================================================================

		# Title of Page 
		label = tk.Label(root, text="Coral PRYSM Model", font=LARGE_FONT)
		label.grid(sticky="E")

		rowIdx = 1

		# =========================================================================================
		# COEFFICIENT ENTRY BOXES
		# =========================================================================================

		# Longitude 
		tk.Label(root, text="Longitude (Range 0 to 360)").grid(row=rowIdx, column=0, sticky="E")
		self.lonEntry = tk.Entry(root)
		self.lonEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Latitude 
		tk.Label(root, text="Latitude (Range -90 to 90)").grid(row=rowIdx, column=0, sticky="E")
		self.latEntry = tk.Entry(root)
		self.latEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Coral Species 
		species = ["Porites_sp", "Porites_lob", "Porites_lut", "Porites_aus", 
		"Montast", "Diploas", "Default"]

		self.v = tk.StringVar()
		self.v.set("default")
		tk.Label(root, text = "Pick a Coral Species:").grid(row=rowIdx, column=0, sticky="E")
		for text in species: 
			b = tk.Radiobutton(root, text=text, variable=self.v, value=text)
			b.grid(row=rowIdx, column=1, sticky="w")
			rowIdx+=1

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

		# Creates example graph (as seen in paper)
		tk.Label(root, text = "Click to see an example graph:").grid(row=rowIdx, column=0, sticky="E")
		exampleGraphButton = tk.Button(root, text="Example Graph", command=self.exampleGraph)
		exampleGraphButton.grid(row=rowIdx, column=1, ipadx=20, ipady=3, sticky="W")
		rowIdx+=1 

		# Allows user to upload data. 
		tk.Label(root, 
			text="Upload a single CSV file with three columns \n and the headers \"TIME\", \"SST\", \"SSS\"."
			).grid(row=rowIdx, columnspan=3, rowspan=3)
		rowIdx += 3
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
		PSGraphButton.grid(row=0, column=12, ipadx=20, ipady=3, sticky="S")

		SensorGraphButton = tk.Button(root, text="Sensor", command=self.generateGraph)
		SensorGraphButton.grid(row=0, column=13, ipadx=20, ipady=3, sticky="S")

		self.f = Figure(figsize=(10,5), dpi=100)
		self.plt = self.f.add_subplot(111)
		self.canvas = FigureCanvasTkAgg(self.f, root)
		self.canvas.get_tk_widget().grid(row=1, column=3, rowspan=16, columnspan=15, sticky="nw")
		self.plt.set_title(r'SENSOR',fontsize=12)
		self.plt.set_xlabel('Time')
		self.plt.set_ylabel('Simulated Coral Data')

		# =========================================================================================
		# SAVE OPTIONS 
		# =========================================================================================	
		saveRowIdx = 18
		tk.Label(root, text = "Save sensor data as:").grid(row=saveRowIdx, column=6, sticky="E")
		dataTxtbutton = tk.Button(root, text=".txt", command=self.saveTxtData)
		dataTxtbutton.grid(row=saveRowIdx, column=8, ipadx=20, ipady=3, sticky="W")
		dataNpybutton = tk.Button(root, text=".npy", command=self.saveNpyData)
		dataNpybutton.grid(row=saveRowIdx, column=9, ipadx=20, ipady=3, sticky="W")
		dataCsvbutton = tk.Button(root, text=".csv", command=self.saveCsvData)
		dataCsvbutton.grid(row=saveRowIdx, column=10, ipadx=20, ipady=3, sticky="W")
		saveRowIdx += 1

		tk.Label(root, text = "Save uncertainty data as:").grid(row=saveRowIdx, column=6, sticky="E")
		errorTxtbutton = tk.Button(root, text=".txt", command=self.saveTxtErrors)
		errorTxtbutton.grid(row=saveRowIdx, column=8, ipadx=20, ipady=3, sticky="W")
		errorCsvbutton = tk.Button(root, text=".csv", command=self.saveCsvErrors)
		errorCsvbutton.grid(row=saveRowIdx, column=9, ipadx=20, ipady=3, sticky="W")
		saveRowIdx += 1

		tk.Label(root, text = "Save power spectrum data as:").grid(row=saveRowIdx, column=6, sticky="E")
		errorTxtbutton = tk.Button(root, text=".txt", command=self.saveTxtPSData)
		errorTxtbutton.grid(row=saveRowIdx, column=8, ipadx=20, ipady=3, sticky="W")
		errorCsvbutton = tk.Button(root, text=".csv", command=self.saveCsvPSData)
		errorCsvbutton.grid(row=saveRowIdx, column=9, ipadx=20, ipady=3, sticky="W")
		saveRowIdx += 1

		tk.Label(root, text = "Save graph as:").grid(row=saveRowIdx, column=6, sticky="E")
		graphPNGbutton = tk.Button(root, text=".png", command=self.savePngGraph)
		graphPNGbutton.grid(row=saveRowIdx, column=8, ipadx=20, ipady=3, sticky="W")
		graphPDFbutton = tk.Button(root, text=".pdf", command=self.savePdfGraph)
		graphPDFbutton.grid(row=saveRowIdx, column=9, ipadx=20, ipady=3, sticky="W")
		saveRowIdx += 1

	"""
	Saves power spectrum data as text file. 
	"""
	def saveTxtPSData(self): 
		if (self.Cd18Of.size != 0):
			file1 = filedialog.asksaveasfilename(initialfile="Cd180f", defaultextension=".txt")
			if file1:
				np.savetxt(file1, self.Cd18Of, newline=" ")
				tk.messagebox.showinfo("Sucess", "Saved power spectrum frequency data")
		if (self.Cd18Opsd_mt.size != 0):
			file2 = filedialog.asksaveasfilename(initialfile="Cd18Opsd_mt", defaultextension=".txt")
			if file2:
				np.savetxt(file2, self.Cd18Opsd_mt, newline=" ")
				tk.messagebox.showinfo("Sucess", "Saved power spectral density data")

	"""
	Saves power spectrum data as csv file. 
	"""
	def saveCsvPSData(self): 
		df = pd.DataFrame({})
		if (self.Cd18Of.size != 0):
			df["Cd180f"] = self.Cd18Of
		if (self.Cd18Opsd_mt.size != 0):
			df["Cd18Opsd_mt"] = self.Cd18Opsd_mt
		file = filedialog.asksaveasfilename(initialfile="PowerSpectrumData.csv", defaultextension=".csv")
		if file:
			df.to_csv(file, index=False)
			tk.messagebox.showinfo("Sucess", "Saved power spectrum data")

	"""
	Saves current graph as pdf.
	"""
	def savePdfGraph(self):
		file = filedialog.asksaveasfilename(initialfile="Figure.pdf", defaultextension=".pdf")
		if file:
			self.f.savefig(file)
			tk.messagebox.showinfo("Sucess", "Saved graph")

	"""
	Saves current graph as png.
	"""
	def savePngGraph(self):
		file = filedialog.asksaveasfilename(initialfile="Figure.png", defaultextension=".png")
		if file:
			self.f.savefig(file)
			tk.messagebox.showinfo("Sucess", "Saved graph")

	"""
	Saves sensor data into a numpy file. 
	"""
	def saveNpyData(self):
		file = filedialog.asksaveasfilename(initialfile="SensorData.npy", defaultextension=".npy")
		if file:
			np.save(file, self.coral)
			tk.messagebox.showinfo("Sucess", "Saved simulated data")

	"""
	Saves sensor data into a text file. 
	"""
	def saveTxtData(self):
		file = filedialog.asksaveasfilename(initialfile="SensorData.txt", defaultextension=".txt")
		if file:
			np.savetxt(file, self.coral, newline=" ")
			tk.messagebox.showinfo("Sucess", "Saved simulated data")
	
	"""
	Saves sensor data into a csv file. 
	"""
	def saveCsvData(self):
		df = pd.DataFrame({"Sensor": self.coral})
		file = filedialog.asksaveasfilename(initialfile="SensorData.csv", defaultextension=".csv")
		if file:
			df.to_csv(file, index=False)
			tk.messagebox.showinfo("Sucess", "Saved uncertainty data")

	"""
	Saves error data into a text file. 
	"""
	def saveCsvErrors(self):
		df = pd.DataFrame({})
		if (self.ageq1.size != 0):
			df["Age_Q1"] = self.ageq1[:,0]
			df["Age_Q2"] = self.ageq1[:,1]

		if (self.gaussq1.size != 0):
			df["GaussianAnalytical_Q1"] = self.gaussq1[:,0]
			df["GaussianAnalytical_Q2"] = self.gaussq1[:,1]
		
		if (self.simpleq1.size != 0 and self.simpleq2.size != 0):
			df["SimpleAnalytical_Q1"] = self.simpleq1
			df["SimpleAnalytical_Q2"] = self.simpleq2
			
		file = filedialog.asksaveasfilename(initialfile="uncertainties.csv", defaultextension=".csv")
		if file:
			df.to_csv(file, index=False)
			tk.messagebox.showinfo("Sucess", "Saved uncertainty data")

	"""
	Saves error data into a text file. 
	"""
	def saveTxtErrors(self):	
		if (self.ageq1.size != 0):
			file1 = filedialog.asksaveasfilename(initialfile="age_uncertainties", defaultextension=".txt")
			if file1:
				np.savetxt(file1, (self.ageq1[:,0], self.ageq1[:,1]), newline=" ")
				tk.messagebox.showinfo("Sucess", 
					"Saved age error data where Q1 is the first row and Q2 is the second row.")
		if (self.gaussq1.size != 0):
			file2 = filedialog.asksaveasfilename(initialfile="gaussian_analytical_uncertainties", defaultextension=".txt")
			if file2:			
				np.savetxt(file2, (self.gaussq1[:,0], self.gaussq1[:,1]), newline=" ")
				tk.messagebox.showinfo("Sucess", 
					"Saved gaussian analytical error data where Q1 is the first row and Q2 is the second row.")
		if (self.simpleq1.size != 0 and self.simpleq2.size != 0):
			file3 = filedialog.asksaveasfilename(initialfile="simple_analytical_uncertainties", defaultextension=".txt")
			if file3:	
				np.savetxt(file3, (self.simpleq1, self.simpleq2), newline=" ")
				tk.messagebox.showinfo("Sucess", 
					"Saved simple analytical error data where Q1 is the first row and Q2 is the second row.")

	"""
	Creates example graph as seen in paper. 
	"""
	def exampleGraph(self): 
		# clear whatever is currently on the canvas 
		self.plt.clear()
		# set newData to true because we are changing the data 
		self.newData = True 
		# set example to true so user will only see example graphs 
		self.example = True 
		# Time values 
		self.time = np.arange(850,1850,1)
		# Get y-values from driver script 
		self.coral = np.load('coral/simulated_coral_d18O.npy')
		self.SST = np.load('coral/test_data_coral/Palmyra_SST_Anomalies_Yearly_850-1850.npy')
		self.SSS = np.load('coral/test_data_coral/Palmyra_SSS_Anomalies_Yearly_850-1850.npy')

		# Reshape coral data for uncertainty calculations 
		self.X = self.coral
		self.X = self.X.reshape(len(self.X),1)

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
			if (self.R1 != float(rateRaw) or self.newData == True):
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
			if (self.R2 != float(sigmaRaw) or self.newData == True):
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
			if (self.R3 != float(sigmaRaw) or self.newData == True):
				sigma = float(sigmaRaw)
				self.R3 = sigma 
				Xn = analytical_error.analytical_error(self.X, sigma, nsamples=100)
				Xn = Xn[:,0,:].reshape(len(self.X), 100) 
				self.gaussq1=mquantiles(Xn,prob=[0.025,0.975],axis=1)
				q2=self.time
				self.plt.fill_between(q2,self.gaussq1[:,0],self.gaussq1[:,1],
					label='100 Gaussian Analytical Error Realizations, CI', facecolor='darkgray',alpha=0.5)

		# Plot the graph
		if (self.PS == True):
			self.ax.remove()
			self.ax2.remove()
			self.ax3.remove()
			self.plt = self.f.add_subplot(111)
			self.canvas = FigureCanvasTkAgg(self.f, root)
			self.canvas.get_tk_widget().grid(row=1, column=3, rowspan=16, columnspan=15, sticky="nw")	
			self.PS = False 

		self.plt.set_title(r'SENSOR')
		self.plt.set_xlabel('Time')
		self.plt.set_ylabel('Simulated Coral Data')		
		self.plt.plot(self.time, self.coral)
		self.canvas.draw()

	"""
	Takes a CSV file and saves the Time, SSS, and SST values. 
	"""
	def uploadData(self):
		# Open the file choosen by the user 
		filename = filedialog.askopenfilename(filetypes = (("csv files","*.csv"),))
		data = np.genfromtxt(filename, delimiter = ",", names=True, dtype=None)
		self.currentFileLabel.configure(text=os.path.basename(filename))
		# Get the entry fields.  
		self.time=data['TIME']
		self.SST=data['SST']
		self.SSS=data['SSS']
		self.newData = True	
		if (self.time.shape != self.SST.shape or self.SST.shape != self.SSS.shape or self.time.shape != self.SSS.shape):
			tk.messagebox.showerror("Error", "Invalid Data: Data inputs are different lengths.")	
			return
	
	"""
	Take the input data and calculate the sensor data. 
	"""	
	def dataPrep(self): 

		# Make sure user has uploaded data 
		if (self.time.size == 0 or self.SST.size == 0 or self.SSS.size == 0):
			tk.messagebox.showerror("Error", "Missing input data")
			return

		# Data has not changed, just return 
		if (self.newData == False): 
			return 

		# Clear whatever is currently on the canvas 
		self.plt.clear()

		# Get the name of the species that was selected. 
		speciesInput = self.v.get()

        # Get the longitude value. 
		lonEntryRaw = self.lonEntry.get()
		if lonEntryRaw == "":
			lon = 197.92
		else:
			lon = float(lonEntryRaw)
			# Longitude must be in the range [0, 360]
			if (lon > 360):
				tk.messagebox.showerror("Error", "Invalid Longitude Value")
				return

		# Get the latitude value. 
		latEntryRaw = self.latEntry.get()
		if latEntryRaw == "": 
			lat = 5.8833
		else:
			lat = float(latEntryRaw)
			# Latitude must be in the range [-90, 90] 
			if (lat < -90):
				tk.messagebox.showerror("Error", "Invalid Latitude Value")
				return

		# Ensure that Sea-Surface Temperature is in degrees C, not Kelvin.
		temp_flag = any(self.SST>200)
		if (temp_flag):
			for i in range(len(self.SST)):
				self.SST[i] = self.SST[i]-274.15

		# Fill coral array with data same size as input vectors.
		self.coral = np.zeros(len(self.time))
		for i in range(len(self.time)):
			self.coral[i] = sensor.pseudocoral(lat,lon,self.SST[i],self.SSS[i],species=speciesInput)

		# Reshape coral data for uncertainty calculations 
		self.X = self.coral
		self.X = self.X.reshape(len(self.X),1)

	"""
	Generates coral data based on input data and model, and graphs result.
	"""
	def generateGraph(self): 

		# user uploads new data or switches from example to non 
		if (self.newData == True or self.example == False):
			self.dataPrep()

		# change from 3 subplots to 1 
		if (self.PS == True): 
			self.ax.remove()
			self.ax2.remove()
			self.ax3.remove()
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
			if (self.R1 != float(rateRaw) or self.newData == True):
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
			if (self.R2 != float(sigmaRaw) or self.newData == True):
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
			if (self.R3 != float(sigmaRaw) or self.newData == True):
				sigma = float(sigmaRaw)
				self.R3 = sigma 
				inputs = len(self.X)
				self.X = self.coral
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
		self.plt.set_ylabel('Simulated Coral Data')
		self.plt.plot(self.time, self.coral)
		self.canvas = FigureCanvasTkAgg(self.f, root)
		self.canvas.get_tk_widget().grid(row=1, column=3, rowspan=16, columnspan=15, sticky="nw")	
		self.canvas.draw()

	"""
	Generates power spectrum of coral model. 
	"""
	def generatePS(self):

		if (self.example == False):
			self.dataPrep()

		if (self.example == True): 
			# Reshape coral data for uncertainty calculations 
			self.X = self.coral
			self.X = self.X.reshape(len(self.X),1)

		# clear whatever is currently on the canvas 
		self.plt.clear()
		self.f.clear()
		self.PS = True 

		if (self.newData == True):
			CST = self.SST - np.mean(self.SST)
			CSS = self.SSS - np.mean(self.SSS)
			Cd18O = self.coral - np.mean(self.coral)

			CSTf, CSTpsd_mt, CSTnu = tsa.multi_taper_psd(CST, Fs=1.0,adaptive=False, jackknife=False)
			CSSf, CSSpsd_mt, CSSnu = tsa.multi_taper_psd(CSS, Fs=1.0,adaptive=False, jackknife=False)
			self.Cd18Of, self.Cd18Opsd_mt, Cd18Onu = tsa.multi_taper_psd(Cd18O, Fs=1.0,adaptive=False, jackknife=False)

		# Get input error rate, if specified 
		rateRaw = self.ageErrorEntry.get()
		if rateRaw == "":
			rate = 0.02
		else:
			rate = float(rateRaw)

		# if rate hasn't changed and already calculated, no need to recalculate 
		if (self.R1 != rate or self.newData == True):
			self.R1 = rate
			if self.Xp.size == 0: 
				tp,self.Xp,tmc = banded.bam_simul_perturb(self.X,self.time,param=[rate,rate],name='poisson',ns=100,resize=0)
			Xpm = self.Xp - np.mean(self.Xp, axis=0)

			# DO SAME CALCULATION IN LOOP FOR ALL AGE UNCERTAINTY VECTORS
			Xpf=np.zeros((len(CSTf),len(Xpm[1])))
			Xpsd_mt=np.zeros((len(CSTf),len(Xpm[1])))
			Xpnu=np.zeros((len(CSTf),len(Xpm[1])))

			# Error here: ValueError: could not broadcast input array from shape (501) into shape (6007)
			for i in range(len(Xpm[1])):
				Xpf[:,i],Xpsd_mt[:,i],Xpnu[:,i]=tsa.multi_taper_psd(Xpm[:,i],Fs=1.0,adaptive=False, jackknife=False)
			# Compute quantiles for spectra
			q1=mquantiles(Xpsd_mt,prob=[0.025,0.975],axis=1)
			# x axis for quantile of spectra:
			q2=Xpf[:,0]			

		# =========================================================================================
		
		self.ax=self.f.add_subplot(311)
		self.ax.spines["top"].set_visible(False)  
		self.ax.spines["right"].set_visible(False) 
		self.ax.loglog(CSTf,CSTpsd_mt, label='SST',color='red')
		self.ax.loglog(CSSf,CSSpsd_mt, label='SSS (PSU)',color='LightSeaGreen')
		self.ax.tick_params(axis="both", which="both", bottom="on", top="off",  
		                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
		self.ax.minorticks_off()
		self.ax.set_title(r'ENVIRONMENT', color='gray')
		self.ax.set_xlabel(r'Period (Years)')
		self.ax.set_ylabel(r'PSD')

		# Change the x-tick labels to something reasonable
		pertick = [500,200,100,50,20,10,8,6,4,2]
		xtick = [500,200,100,50,20,10,8,6,4,2]
		for i in range(len(pertick)):
			xtick[i] = 1.0/pertick[i]
		pertick_labels=(['500', '200', '100' ,'50', '20','10','8' ,'6','4', '2'])

		self.ax.set_xticks(xtick)
		self.ax.set_xticklabels(pertick_labels)
		self.ax.grid('on',axis='y',color='DimGray')
		self.ax.set_xlim([1./550.,0.4])
		self.ax.set_ylim([1e-3, 1e1])
		self.ax.legend(loc=3,fontsize=11,frameon=False)

		# =========================================================================================

		self.ax2=self.f.add_subplot(312)
		self.ax2.spines["top"].set_visible(False)  
		self.ax2.spines["right"].set_visible(False)  
		self.ax2.loglog(self.Cd18Of,self.Cd18Opsd_mt, label=r'Coral $\delta^{18}O_{C}$',color='DarkOrange')
		self.ax2.set_title(r'SENSOR', color='gray')
		self.ax2.set_xlabel(r'Period (Years)')
		self.ax2.set_ylabel(r'PSD')
		self.ax2.tick_params(axis="both", which="both", bottom="on", top="off",  
		                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
		self.ax2.minorticks_off()

		self.ax2.grid('on',axis='y',color='DimGray')
		self.ax2.set_xticks(xtick)
		self.ax2.set_xticklabels(pertick_labels)
		self.ax2.set_xlim([1./550.,0.4])
		self.ax2.set_ylim([1e-3, 1])
		self.ax2.legend(loc=3,fontsize=11,frameon=False)

		# =========================================================================================

		self.ax3=self.f.add_subplot(313)
		self.ax3.spines["top"].set_visible(False)  
		self.ax3.spines["right"].set_visible(False)  
		# Observation purturbed age ensemble here.
		self.ax3.fill_between(q2,q1[:,0],q1[:,1],label='1000 Age-Perturbed Realizations (4%), CI',facecolor='gray',alpha=0.5)
		self.ax3.loglog(self.Cd18Of,self.Cd18Opsd_mt, label=r'Coral $\delta^{18}O_{C}$',color='DarkOrange')
		self.ax3.tick_params(axis="both", which="both", bottom="on", top="off",  
		                labelbottom="on", left="on", right="off", labelleft="on",direction="out")  
		self.ax3.minorticks_off()
		self.ax3.set_title(r'OBSERVATION', color='gray')
		self.ax3.set_xlabel(r'Period (Years)')
		self.ax3.set_ylabel(r'PSD')

		self.ax3.set_xticks(xtick)
		self.ax3.set_xticklabels(pertick_labels)
		self.ax3.grid('on',axis='y',color='DimGray')
		self.ax3.set_xlim([1./550.,0.4])
		self.ax3.set_ylim([1e-3, 1])
		self.ax3.legend(loc=3,fontsize=11,frameon=False)
	
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
			self.ax3.cla()
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
		self.lonEntry.delete(0, tk.END)
		self.latEntry.delete(0, tk.END)		
		self.v.set(None)


"""
Initialize Application.
"""
if __name__ == "__main__":
    root = tk.Tk()
    Grapher(root)
    root.geometry("1600x700+0+50")
    root.mainloop()

