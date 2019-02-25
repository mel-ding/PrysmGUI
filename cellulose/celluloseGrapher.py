import sensor as sensor

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
		self.parent.title("Cellulose Grapher")

		# Title of Page
		label = tk.Label(root, text="Cellulose Graph Page", font=LARGE_FONT)
		label.grid(padx=(25,0))

		rowIdx = 1

		# COEFFICIENT ENTRY BOXES =================================================================

		# Time (vector array, Years)
		tk.Label(root, text = "Time (Years)").grid(row=1, column=0)
		self.lonEntry = tk.Entry(root)
		self.lonEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Temperature (Kelvin)
		tk.Label(root, text = "Temperature").grid(row=2, column=0)
		self.latEntry = tk.Entry(root)
		self.latEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Precipitation mm/day
		tk.Label(root, text = "Precipitation (mm/day)").grid(row=3, column=0)
		self.latEntry = tk.Entry(root)
		self.latEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# Relative Humidity %
		tk.Label(root, text = "Relative Humidity (%)").grid(row=4, column=0)
		self.latEntry = tk.Entry(root)
		self.latEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# d180s isotope ratio of soil water
		tk.Label(root, text = "d180s (Isotope ratio of Soil Water)").grid(row=5, column=0)
		self.latEntry = tk.Entry(root)
		self.latEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# d180p isotope ratio of precipitation
		tk.Label(root, text = "d180p (Isotope ratio of Precipitation)").grid(row=6, column=0)
		self.latEntry = tk.Entry(root)
		self.latEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		# d180v isotope ratio of Ambient vapor at surface layer
		tk.Label(root, text = "d180v (Isotope ratio of Ambient Vapor)").grid(row=7, column=0)
		self.latEntry = tk.Entry(root)
		self.latEntry.grid(row=rowIdx, column=1)
		rowIdx += 1

		Models = [("Roden", 0), ("Evans (Default)", 1)]

		self.v = tk.StringVar()
		self.v.set("default")
		tk.Label(root, text = "Pick a Cellulose Model:").grid(row=8, column=0)
		for text, mode in Models:
			b = tk.Radiobutton(root, text=text, variable=self.v, value=mode)
			b.grid(row=rowIdx, column=1, sticky="w")
			rowIdx+=1


		# GRAPHER AND BUTTONS ====================================================================

		# Creates example graph (as seen in paper)
		tk.Label(root, text = "Click to see an example graph:").grid(row=rowIdx, column=0)
		exampleGraphButton = tk.Button(root, text="Example Graph", command=self.exampleGraph)
		exampleGraphButton.grid(row=rowIdx, column=1)
		rowIdx+=1

		# Creates graph with given inputs
		tk.Label(root,
			text="Upload a single CSV file with three columns \n and the headers \"TIME\", \"SST\", \"SSS\"."
			).grid(row=rowIdx, columnspan=3, rowspan=3)
		rowIdx += 3
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
		yVals = np.load('simulated_coral_d18O.npy')
		# ASK: WHAT IS CORAL AGE PERTURBED? error bars
		errorVals = np.load('coral_age_perturbed.npy')

		# Plot the graph
		self.a.plot(xVals, yVals)
		self.a.set_xlabel('Time')
		self.a.set_ylabel('Simulated Coral Data')
		self.canvas.draw()

	"""
	Takes a CSV file with Time, SSS, and SST values and generates
	coral data based on the model.
	"""
	def generateGraph(self):
		# Open the file choosen by the user
		filename = filedialog.askopenfilename(filetypes = (("csv files","*.csv"),))
		data = np.genfromtxt(filename, delimiter = ",", names=True, dtype=None)
		# Get the entry fields.
		time=data['TIME']
		SST=data['SST']
		SSS=data['SSS']

		# Clear whatever is currently on the canvas
		self.a.clear()

		speciesInput = self.v.get()
		print(speciesInput)

        # Longitude   [lon]       [0, 360]
		lonEntryRaw = self.lonEntry.get()
		if lonEntryRaw == "":
			lon = 197.92
		else:
			lon = float(lonEntryRaw)
			if (lon<0.):
				lon = lon+360.

		# Latitude    [lat]       [-90, 90]
		latEntryRaw = self.latEntry.get()
		if latEntryRaw == "":
			lat = 5.8833
		else:
			lat = float(latEntryRaw)
			if (lat>90.):
				lat = lat-90.

		# Ensure that Sea-Surface Temperature is in degrees C, not Kelvin.
		temp_flag = any(SST>200)
		if (temp_flag):
			for i in range(len(SST)):
				SST[i] = SST[i]-274.15

		# Fill coral array with data same size as input vectors.
		coral = np.zeros(len(time))
		for i in range(len(time)):
			coral[i] = sensor.pseudocoral(lat,lon,SST[i],SSS[i],species=speciesInput)

		# Plot the graph
		self.a.plot(time, coral)
		self.a.set_xlabel('Time')
		self.a.set_ylabel('Simulated Coral Data')
		self.canvas.draw()

	"""
	Clears the content in the graph.
	TODO: Make clear button clear entry field contents
	"""
	def clear(self):
		self.a.clear()
		self.canvas.draw()

"""
Initialize Application.
"""
if __name__ == "__main__":
    root = tk.Tk()
    Grapher(root)
    root.geometry("1375x550+50+100")
    root.mainloop()
