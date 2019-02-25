import os
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

		models = ["Roden", "Evans"]

		self.v = tk.StringVar()
		self.v.set("default")
		tk.Label(root, text = "Pick a Cellulose Model:").grid(row=3, column=0)
		for text in models:
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
			text="Upload a single CSV file with seven columns and the headers \n \"TIME\", \"TEMP\", \"PRECIP\", \"RELHUM\", \"D180S\", \"D180P\", \"D180V\".").grid(row=rowIdx, columnspan=3, rowspan=3)
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
		time = np.arange(1000,2005,1)
		# Get y-values from driver script
		yVals = np.load('../dcell_Xn.npy')
		# ASK: WHAT IS CELLULOSE AGE PERTURBED? error bars
		# errorVals = np.load('coral/coral_age_perturbed.npy') # To be changed

		# Plot the graph
		self.a.plot(time, yVals)
		self.a.set_xlabel('Time')
		self.a.set_ylabel('Simulated Cellulose Data')
		self.canvas.draw()

	"""
	Takes a CSV file and saves all arrays needed for Cellulose models.
	"""
	def uploadData(self):
		# Open the file choosen by the user
		filename = filedialog.askopenfilename(filetypes = (("csv files","*.csv"),))
		data = np.genfromtxt(filename, delimiter = ",", names=True, dtype=None)
		self.currentFileLabel.configure(text=os.path.basename(filename))

		# Get the entry fields.
		self.time=data['TIME']
		self.temp=data['TEMP']
		self.precip=data['PRECIP']
		self.rh=data['RELHUM']
		self.d180s=data['D180S']
		self.d180p=data['D180P']
		self.d180v=data['D180V']

	"""
	Generates coral data based on input data and model, and graphs result.
	"""
	def generateGraph(self):
		# Clear whatever is currently on the canvas
		self.a.clear()

		# Get the name of the model that was selected.
		model = self.v.get()
		flag = 0 if model == "Roden" else 1

		# Fill coral array with data same size as input vectors.
		cellulose = np.zeros(len(self.time))
		for i in range(len(self.time)):
			cellulose[i] = sensor.cellulose_sensor(time[i], self.temp[i], self.precip[i],
			self.rh[i], self.d180p[i], self.d180s[i], self.d180v, flag)[0]

		# Plot the graph
		self.a.plot(self.time, cellulose)
		self.a.set_xlabel('Time')
		self.a.set_ylabel('Simulated Cellulose Data')
		self.canvas.draw()

	"""
	Clears the content in the graph.
	TODO: Make clear button clear entry field contents
	"""
	def clear(self):
		# self.lonEntry.delete(0, tk.END)
		# self.latEntry.delete(0, tk.END)
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
