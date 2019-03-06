import os
import psm.coral.sensor as sensor
import psm.agemodels.banded as banded
import psm.aux_functions.analytical_err_simple as analytical_err_simple
import psm.aux_functions.analytical_error as analytical_error
from scipy.stats.mstats import mquantiles

import tkinter as tk
from tkinter import ttk, filedialog
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

LARGE_FONT = ("Verdana", 25) 


class Grapher(tk.Frame):

	def __init__(self, parent, *args, **kwargs):
		tk.Frame.__init__(self, parent, *args, **kwargs)
		self.parent = parent
		self.parent.title("Coral Grapher")

		# Title of Page 
		label = tk.Label(root, text="Coral Graph Page", font=LARGE_FONT)
		label.grid(padx=(25,0))

		rowIdx = 1

		# COEFFICIENT ENTRY BOXES =================================================================

		# Longitude 
		tk.Label(root, text="Longitude [0, 360]").grid(row=rowIdx, column=0)
		self.lonEntry = tk.Entry(root)
		self.lonEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Latitude 
		tk.Label(root, text="Latitude [-90, 90]").grid(row=rowIdx, column=0)
		self.latEntry = tk.Entry(root)
		self.latEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Coral Species 
		species = ["Porites_sp", "Porites_lob", "Porites_lut", "Porites_aus", 
		"Montast", "Diploas", "Default"]

		self.v = tk.StringVar()
		self.v.set("default")
		tk.Label(root, text = "Pick a Coral Species:").grid(row=rowIdx, column=0)
		for text in species: 
			b = tk.Radiobutton(root, text=text, variable=self.v, value=text)
			b.grid(row=rowIdx, column=1, sticky="w")
			rowIdx+=1

		# ERROR OPTIONS ==========================================================================

		# TODO: ask if want checkboxes or radio buttons? 
		# is showing it on the right okay?

		# Age Uncertainties
		self.ageErrorVal = tk.BooleanVar() 
		self.ageErrorVar = tk.Checkbutton(root, variable=self.ageErrorVal)
		self.ageErrorVar.grid(row=rowIdx, column=2)
		tk.Label(root, text = "Add age uncertainties with error rate of:").grid(row=rowIdx, column=0)
		self.ageErrorEntry = tk.Entry(root)
		self.ageErrorEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Analytical Errors
		self.altErrorVal = tk.BooleanVar() 
		self.altErrorVar = tk.Checkbutton(root, variable=self.altErrorVal)
		self.altErrorVar.grid(row=rowIdx, column=2)
		tk.Label(root, text = "Add simple analytical error with precision of:").grid(row=rowIdx, column=0)
		self.altErrorEntry = tk.Entry(root)
		self.altErrorEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# GRAPHER AND BUTTONS ====================================================================

		# Creates example graph (as seen in paper)
		tk.Label(root, text = "Click to see an example graph:").grid(row=rowIdx, column=0)
		exampleGraphButton = tk.Button(root, text="Example Graph", command=self.exampleGraph)
		exampleGraphButton.grid(row=rowIdx, column=1)
		rowIdx+=1 

		# Allows user to upload data. 
		tk.Label(root, 
			text="Upload a single CSV file with three columns \n and the headers \"TIME\", \"SST\", \"SSS\"."
			).grid(row=rowIdx, columnspan=3, rowspan=3)
		rowIdx += 3
		tk.Label(root, text = "Click to upload your data:").grid(row=rowIdx, column=0)
		graphButton = tk.Button(root, text="Upload Data", command=self.uploadData)
		graphButton.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Shows the name of the current uploaded file, if any. 
		tk.Label(root, text="Current File Uploaded:").grid(row=rowIdx, column=0)
		self.currentFileLabel = tk.Label(root, text="No file")
		self.currentFileLabel.grid(row=rowIdx, column=1, columnspan=2)
		rowIdx+=1

		# Creates graph with given inputs
		tk.Label(root, text = "Click to create your own graph:").grid(row=rowIdx, column=0)
		graphButton = tk.Button(root, text="Generate Graph", command=self.generateGraph)
		graphButton.grid(row=rowIdx, column=1)
		rowIdx+=1

		# Clears graph
		tk.Label(root, text = "Click to clear the graph:").grid(row=rowIdx, column=0)
		clearButton = tk.Button(root, text="Clear", command=self.clear)
		clearButton.grid(row=rowIdx, column=1)

		f = Figure(figsize=(10,5), dpi=100)
		self.plt = f.add_subplot(111)
		self.canvas = FigureCanvasTkAgg(f, root)
		self.canvas.get_tk_widget().grid(row=0, column=3, rowspan=16, columnspan=15, sticky="nw")

	"""
	Creates example graph as seen in paper. 
	"""
	def exampleGraph(self): 
		# clear whatever is currently on the canvas 
		self.plt.clear()
		# Time values 
		self.time = np.arange(850,1850,1)
		# Get y-values from driver script 
		coral = np.load('coral/simulated_coral_d18O.npy')
		# Plot age uncertainties - if selected 
		if (self.ageErrorVal.get()):
			# Get the input error rate value
			rateRaw = self.ageErrorEntry.get()
			if rateRaw == "":
				tk.messagebox.showerror("Error", "Invalid Rate Value")
				return
			rate = float(rateRaw)
			# Calculate the age uncertanties
			X = coral
			X = X.reshape(len(X),1)
			tp, Xp, tmc = banded.bam_simul_perturb(X, self.time, param=[rate, rate])
			q1=mquantiles(Xp,prob=[0.025,0.975],axis=1)
			q2=self.time
			self.plt.fill_between(q2,q1[:,0],q1[:,1],label='1000 Age-Perturbed Realizations, CI',facecolor='gray',alpha=0.5)

		# Plot the graph
		self.plt.plot(self.time, coral)
		self.plt.set_xlabel('Time')
		self.plt.set_ylabel('Simulated Coral Data')
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
		# TODO: Check that this error catching actually works
		if (self.time.shape != self.SST.shape or self.SST.shape != self.SSS.shape or self.time.shape != self.SSS.shape):
			tk.messagebox.showerror("Error", "Invalid Data: Data inputs are different lengths.")
		# print(self.time.shape, self.SST.shape, self.SSS.shape)	
	
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
		coral = np.zeros(len(self.time))
		for i in range(len(self.time)):
			coral[i] = sensor.pseudocoral(lat,lon,self.SST[i],self.SSS[i],species=speciesInput)

		# Plot age uncertainties - if selected 
		if (self.ageErrorVal.get()):
			# Get the input error rate value
			rateRaw = self.ageErrorEntry.get()
			if rateRaw == "":
				tk.messagebox.showerror("Error", "Invalid Rate Value")
				return
			rate = float(rateRaw)
			# Calculate the age uncertanties
			X = coral
			X = X.reshape(len(X),1)
			tp, Xp, tmc = banded.bam_simul_perturb(X, self.time, param=[rate, rate])
			q1=mquantiles(Xp,prob=[0.025,0.975],axis=1)
			q2=self.time
			self.plt.fill_between(q2,q1[:,0],q1[:,1],label='1000 Age-Perturbed Realizations, CI',facecolor='gray',alpha=0.5)

		# Plot analytical errors - if selected
		if (self.altErrorVal.get()):
			# Get the input error rate value
			sigmaRaw = self.altErrorEntry.get()
			if sigmaRaw == "":
				tk.messagebox.showerror("Error", "Invalid Rate Value")
				return
			sigma = float(sigmaRaw)
			X = coral
			X = X.reshape(len(X),1)
			q1, q2 = analytical_err_simple.analytical_err_simple(X,sigma)
			q1 = q1.reshape(len(self.time))
			q2 = q2.reshape(len(self.time))
			self.plt.fill_between(self.time,q1,q2,label='100 Analytical Error Realizations, CI',facecolor='darkgray',alpha=0.5)

		"""
		for analytical errors need to use quantile
		q1=mquantiles(Xp,prob=[0.025,0.975],axis=1)
		q2=self.time

		add error rates to clear function 

		download options - data and figure 
		"""

		# Plot the graph
		self.plt.plot(self.time, coral)
		self.plt.set_xlabel('Time')
		self.plt.set_ylabel('Simulated Coral Data')
		self.canvas.draw()

	"""
	Clears the content in the graph. 
	TODO: Make clear button clear entry field contents
	"""
	def clear(self):
		self.lonEntry.delete(0, tk.END)
		self.latEntry.delete(0, tk.END)
		self.v.set(None)
		self.plt.clear()
		self.canvas.draw()

"""
Initialize Application.
"""
if __name__ == "__main__":
    root = tk.Tk()
    Grapher(root)
    root.geometry("1440x620+5+100")
    root.mainloop()
