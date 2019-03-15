import os
import psm.coral.sensor as sensor
import psm.agemodels.banded as banded
import psm.aux_functions.analytical_err_simple as analytical_err_simple
import psm.aux_functions.analytical_error as analytical_error
from scipy.stats.mstats import mquantiles

import tkinter as tk
from tkinter import ttk, filedialog
import numpy as np
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
		self.parent.title("Coral PRSYM Models")

		# Initialize empty arrays for data saving 
		self.ageq1 = np.array([])
		self.gaussq1 = np.array([])
		self.simpleq1 = np.array([])
		self.simpleq2 = np.array([])

		# =========================================================================================
		# SENSOR DATA
		# =========================================================================================

		# Title of Page 
		label = tk.Label(root, text="Coral Sensor Data", font=LARGE_FONT)
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
		# SAVE OPTIONS 
		# =========================================================================================	
		tk.Label(root, text = "Save sensor data as:").grid(row=18, column=6, sticky="E")
		dataTxtbutton = tk.Button(root, text=".txt", command=self.saveTxtData)
		dataTxtbutton.grid(row=18, column=8, ipadx=20, ipady=3, sticky="W")
		dataNpybutton = tk.Button(root, text=".npy", command=self.saveNpyData)
		dataNpybutton.grid(row=18, column=9, ipadx=20, ipady=3, sticky="W")
		dataCsvbutton = tk.Button(root, text=".npy", command=self.saveCsvData)
		dataCsvbutton.grid(row=18, column=10, ipadx=20, ipady=3, sticky="W")

		tk.Label(root, text = "Save error data as:").grid(row=19, column=6, sticky="E")
		errorTxtbutton = tk.Button(root, text=".txt", command=self.saveTxtErrors)
		errorTxtbutton.grid(row=19, column=8, ipadx=20, ipady=3, sticky="W")
		errorCsvbutton = tk.Button(root, text=".csv", command=self.saveCsvErrors)
		errorCsvbutton.grid(row=19, column=9, ipadx=20, ipady=3, sticky="W")

		tk.Label(root, text = "Save graph as:").grid(row=20, column=6, sticky="E")
		graphPNGbutton = tk.Button(root, text=".png", command=self.savePngGraph)
		graphPNGbutton.grid(row=20, column=8, ipadx=20, ipady=3, sticky="W")
		graphPDFbutton = tk.Button(root, text=".pdf", command=self.savePdfGraph)
		graphPDFbutton.grid(row=20, column=9, ipadx=20, ipady=3, sticky="W")

		# =========================================================================================
		# GRAPH
		# =========================================================================================
		self.f = Figure(figsize=(10,5), dpi=100)
		self.plt = self.f.add_subplot(111)
		self.canvas = FigureCanvasTkAgg(self.f, root)
		self.canvas.get_tk_widget().grid(row=0, column=3, rowspan=16, columnspan=15, sticky="nw")
		self.plt.set_xlabel('Time')
		self.plt.set_ylabel('Simulated Coral Data')

	"""
	Saves current graph as pdf.
	"""
	def savePdfGraph(self):
		self.f.savefig('graph.pdf')
		tk.messagebox.showinfo("Sucess", "Saved graph as graph.pdf")

	"""
	Saves current graph as png.
	"""
	def savePngGraph(self):
		self.f.savefig('graph.png')
		tk.messagebox.showinfo("Sucess", "Saved graph as graph.png")

	"""
	Saves sensor data into a numpy file. 
	"""
	def saveNpyData(self):
		np.save("simulated_coral_d18O.npy", self.coral)
		tk.messagebox.showinfo("Sucess", "Saved simulated data as 'simulated_coral_d18O.npy'")

	"""
	Saves sensor data into a text file. 
	"""
	def saveTxtData(self):
		np.savetxt("simulated_coral_d18O.txt", self.coral, newline=" ")
		tk.messagebox.showinfo("Sucess", "Saved simulated data as 'simulated_coral_d18O.txt'")

	"""
	Saves sensor data into a text file. 
	"""
	def saveTxtData(self):
		np.savetxt("simulated_coral_d18O.csv", self.coral, newline=" ")
		tk.messagebox.showinfo("Sucess", "Saved simulated data as 'simulated_coral_d18O.csv'")
	
	"""
	Saves sensor data into a csv file. 
	"""
	def saveCsvData(self):
		df = pd.DataFrame({"Sensor": self.coral})
		df.to_csv("sensor_data.csv", index=False)
		tk.messagebox.showinfo("Sucess", "Saved error data as 'sensor_data.csv'")


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
			
		df.to_csv("error_data.csv", index=False)
		tk.messagebox.showinfo("Sucess", "Saved error data as 'error_data.csv'")

	"""
	Saves error data into a text file. 
	"""
	def saveTxtErrors(self):	
		if (self.ageq1.size != 0):
			np.savetxt("age_errors.txt", (self.ageq1[:,0], self.ageq1[:,1]), newline=" ")
			tk.messagebox.showinfo("Sucess", 
				"Saved age error data as 'age_errors.txt' where Q1 is the first row and Q2 is the second row.")
		if (self.gaussq1.size != 0):
			np.savetxt("gaussian_errors.txt", (self.gaussq1[:,0], self.gaussq1[:,1]), newline=" ")
			tk.messagebox.showinfo("Sucess", 
				"Saved gaussian analytical error data as 'gaussian_errors.txt' where Q1 is the first row and Q2 is the second row.")
		if (self.simpleq1.size != 0 and self.simpleq2.size != 0):
			np.savetxt("simple_errors.txt", (self.simpleq1, self.simpleq2), newline=" ")
			tk.messagebox.showinfo("Sucess", 
				"Saved simple analytical error data as 'simple_errors.txt' where Q1 is the first row and Q2 is the second row.")

	"""
	Creates example graph as seen in paper. 
	"""
	def exampleGraph(self): 
		# clear whatever is currently on the canvas 
		self.plt.clear()
		# Time values 
		self.time = np.arange(850,1850,1)
		# Get y-values from driver script 
		self.coral = np.load('coral/simulated_coral_d18O.npy')

		# Reshape coral data for uncertainty calculations 
		X = self.coral
		X = X.reshape(len(X),1)

		# Plot age uncertainties - if selected 
		if (self.ageErrorVal.get()):
			# Get input error rate
			rateRaw = self.ageErrorEntry.get()
			if rateRaw == "":
				tk.messagebox.showerror("Error", "Invalid Rate Value")
				return
			rate = float(rateRaw)
			# Calculate the age uncertanties
			tp, Xp, tmc = banded.bam_simul_perturb(X, self.time, param=[rate, rate])
			self.ageq1=mquantiles(Xp,prob=[0.025,0.975],axis=1)
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
			sigma = float(sigmaRaw)
			self.simpleq1, self.simpleq2 = analytical_err_simple.analytical_err_simple(X,sigma)
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
			sigma = float(sigmaRaw)
			inputs = len(X)
			X = self.coral
			X = X.reshape(inputs,1)
			## TO ASK: HAD TO CHANGE NSAMPLES TO MATCH NUMBER OF CORAL DATA POINTS
			Xn = analytical_error.analytical_error(X, sigma, inputs)
			Xn = Xn[:,0,:].reshape(inputs, inputs) 
			self.gaussq1=mquantiles(Xn,prob=[0.025,0.975],axis=1)
			q2=self.time
			self.plt.fill_between(q2,self.gaussq1[:,0],self.gaussq1[:,1],
				label='100 Gaussian Analytical Error Realizations, CI', facecolor='darkgray',alpha=0.5)

		# Plot the graph
		self.plt.plot(self.time, self.coral)
		self.canvas.draw()
		self.plt.set_xlabel('Time')
		self.plt.set_ylabel('Simulated Coral Data')

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
		if (self.time.shape != self.SST.shape or self.SST.shape != self.SSS.shape or self.time.shape != self.SSS.shape):
			tk.messagebox.showerror("Error", "Invalid Data: Data inputs are different lengths.")	
	
	"""
	Generates coral data based on input data and model, and graphs result.
	"""
	def generateGraph(self): 
		# Make sure user has uploaded data 
		if (self.time.size == 0 or self.SST.size == 0 or self.SSS.size == 0):
			tk.messagebox.showerror("Error", "Missing input data")
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
		X = self.coral
		X = X.reshape(len(X),1)

		# Plot age uncertainties - if selected 
		if (self.ageErrorVal.get()):
			# Get input error rate
			rateRaw = self.ageErrorEntry.get()
			if rateRaw == "":
				tk.messagebox.showerror("Error", "Invalid Rate Value")
				return
			rate = float(rateRaw)
			# Calculate the age uncertanties
			tp, Xp, tmc = banded.bam_simul_perturb(X, self.time, param=[rate, rate])
			self.ageq1=mquantiles(Xp,prob=[0.025,0.975],axis=1)
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
			sigma = float(sigmaRaw)
			self.simpleq1, self.simpleq2 = analytical_err_simple.analytical_err_simple(X,sigma)
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
			sigma = float(sigmaRaw)
			inputs = len(X)
			X = self.coral
			X = X.reshape(inputs,1)
			# print(X.shape) # (1001, 1)
			## TO ASK: HAD TO CHANGE NSAMPLES TO MATCH NUMBER OF CORAL DATA POINTS
			Xn = analytical_error.analytical_error(X, sigma, inputs)
			# print(Xn.shape) # (1001, 1001, 1000)
			Xn = Xn[:,0,:].reshape(inputs, inputs) 
			self.gaussq1=mquantiles(Xn,prob=[0.025,0.975],axis=1)
			q2=self.time
			self.plt.fill_between(q2,self.gaussq1[:,0],self.gaussq1[:,1],
				label='100 Gaussian Analytical Error Realizations, CI', facecolor='darkgray',alpha=0.5)

		# Plot the graph
		self.plt.plot(self.time, self.coral)
		self.canvas.draw()
		self.plt.set_xlabel('Time')
		self.plt.set_ylabel('Simulated Coral Data')

	"""
	Clears all content in the graph. 
	"""
	def clearGraph(self):
		self.plt.clear()
		self.canvas.draw()

	"""
	Clears entries the user has input. 
	"""
	def clearEntries(self):
		self.ageErrorVal.set(0)
		self.simAltErrorVal.set(0)
		self.altErrorVal.set(0)
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
    root.geometry("1440x650+5+100")
    root.mainloop()

