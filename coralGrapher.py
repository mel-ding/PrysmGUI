import os
import psm.coral.sensor as sensor

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
		tk.Label(root, text = "Longitude").grid(row=1, column=0)
		self.lonEntry = tk.Entry(root)
		self.lonEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Latitude 
		tk.Label(root, text = "Latitude").grid(row=2, column=0)
		self.latEntry = tk.Entry(root)
		self.latEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		species = ["Porites_sp", "Porites_lob", "Porites_lut", "Porites_aus", 
		"Montast", "Diploas", "Default"]

		self.v = tk.StringVar()
		self.v.set("default")
		tk.Label(root, text = "Pick a Coral Species:").grid(row=3, column=0)
		for text in species: 
			b = tk.Radiobutton(root, text=text, variable=self.v, value=text)
			b.grid(row=rowIdx, column=1, sticky="w")
			rowIdx+=1


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
		self.currentFileLabel.grid(row=rowIdx, column=1)
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
		self.a = f.add_subplot(111)
		self.canvas = FigureCanvasTkAgg(f, root)
		self.canvas.get_tk_widget().grid(row=0, column=3, rowspan=16, columnspan=10, sticky="nw")

	"""
	Creates example graph as seen in paper. 
	"""
	def exampleGraph(self): 
		# clear whatever is currently on the canvas 
		self.a.clear()
		# Time values 
		xVals = np.arange(850,1850,1)
		# Get y-values from driver script 
		yVals = np.load('coral/simulated_coral_d18O.npy')
		# ASK: WHAT IS CORAL AGE PERTURBED? error bars 
		errorVals = np.load('coral/coral_age_perturbed.npy')

		# Plot the graph
		self.a.plot(xVals, yVals)
		self.a.set_xlabel('Time')
		self.a.set_ylabel('Simulated Coral Data')
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

	"""
	Generates coral data based on input data and model, and graphs result.
	"""
	def generateGraph(self): 
		# Clear whatever is currently on the canvas 
		self.a.clear()

		# Get the name of the species that was selected. 
		speciesInput = self.v.get()

        # Get the longitude value. 
		lonEntryRaw = self.lonEntry.get()
		if lonEntryRaw == "":
			lon = 197.92
		else:
			lon = float(lonEntryRaw)
			# Longitude must be in the range [0, 360] or [-180, 180]
			if (lon > 360 or lon < -180):
				tk.messagebox.showerror("Error", "Invalid Longitude Value")
				return
			if (lon < 0.):
				lon += 360.

		# Get the latitude value. 
		latEntryRaw = self.latEntry.get()
		if latEntryRaw == "": 
			lat = 5.8833
		else:
			lat = float(latEntryRaw)
			# Latitude must be in the range [-90, 90] or [0, 180]
			if (lat < -90 or lat > 180):
				tk.messagebox.showerror("Error", "Invalid Latitude Value")
				return
			if (lat > 90.):
				lat -= 90.

		# Ensure that Sea-Surface Temperature is in degrees C, not Kelvin.
		temp_flag = any(self.SST>200)
		if (temp_flag):
			for i in range(len(self.SST)):
				self.SST[i] = self.SST[i]-274.15

		# Fill coral array with data same size as input vectors.
		coral = np.zeros(len(self.time))
		for i in range(len(self.time)):
			coral[i] = sensor.pseudocoral(lat,lon,self.SST[i],self.SSS[i],species=speciesInput)

		# Plot the graph
		self.a.plot(self.time, coral)
		self.a.set_xlabel('Time')
		self.a.set_ylabel('Simulated Coral Data')
		self.canvas.draw()

	"""
	Clears the content in the graph. 
	TODO: Make clear button clear entry field contents
	"""
	def clear(self):
		self.lonEntry.delete(0, tk.END)
		self.latEntry.delete(0, tk.END)
		self.v.set(None)
		self.a.clear()
		self.canvas.draw()

"""
Initialize Application.
"""
if __name__ == "__main__":
    root = tk.Tk()
    Grapher(root)
    root.geometry("1400x590+30+100")
    root.mainloop()
