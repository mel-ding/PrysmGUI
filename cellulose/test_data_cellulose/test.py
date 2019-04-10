# from tkinter import *

# root = Tk()
# root.geometry("600x400")

# left = Frame(root, borderwidth=2, relief="solid")
# right = Frame(root, borderwidth=2, relief="solid")
# container = Frame(left, borderwidth=2, relief="solid")
# box1 = Frame(right, borderwidth=2, relief="solid")
# box2 = Frame(right, borderwidth=2, relief="solid")

# label1 = Label(container, text="I could be a canvas, but I'm a label right now")
# label2 = Label(left, text="I could be a button")
# label3 = Label(left, text="So could I")
# label4 = Label(box1, text="I could be your image")
# label5 = Label(box2, text="I could be your setup window")

# left.pack(side="left", expand=True, fill="both")
# right.pack(side="right", expand=True, fill="both")
# container.pack(expand=True, fill="both", padx=5, pady=5)
# box1.pack(expand=True, fill="both", padx=10, pady=10)
# box2.pack(expand=True, fill="both", padx=10, pady=10)

# label1.pack()
# label2.pack()
# label3.pack()
# label4.pack()
# label5.pack()

# root.mainloop()

# from tkinter import *

# m = PanedWindow(orient=VERTICAL)
# m.pack(fill=BOTH, expand=1)

# top = Label(m, text="top pane")
# m.add(top)

# bottom = Label(m, text="bottom pane")
# m.add(bottom)

# mainloop()

#################################################################################### 

import numpy as np 
import pandas as pd
import os

cwd = os.getcwd()

# for file in os.listdir("cellulose/test_data_cellulose"):
#     print(file)

print("!!!!!!!!!!") 

soil = np.load(cwd + "/d180_soil.npy")
precip = np.load("precipitation.npy")
hum = np.load("relative_hum.npy")
vapor = np.load("d180_vapor.npy")
temp = np.load("temperature.npy")
d180precip = np.load("d180_precip.npy")

timeVals = np.linspace(850, 1850, 1001)
df = pd.DataFrame({"TIME": timeVals, "SOIL":soil, "PRECIPTATION":precip, 
	"HUMIDITY":relative_hum, "VAPOR": vapor, "TEMPERATURE": temp, 
	"D180-PRECIPITATION":d180precip})

df.to_csv("test_data_cellulose.csv", index=False)

#################################################################################### 

# import sys
# import matplotlib
# import matplotlib.pyplot as plt
# matplotlib.use("TkAgg")
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# from matplotlib.figure import Figure
# # Tkinter is for python 2; tkinter is for python 3
# if sys.version_info[0] < 3:
#     import Tkinter as tk
#     import tkMessageBox, tkFileDialog

# else:
#     import tkinter as tk
#     from tkinter import messagebox as tkMessageBox
#     from tkinter import filedialog as tkFileDialog

# class MainApp(tk.Frame):

#     def __init__(self, parent):
#         tk.Frame.__init__(self, parent)
#         self.parent = parent
#         self.parent.title('App')
#         # call the widgets
#         self.okButton()
#         self.quitButton()
#         self.readDataButton()
#         self.clearDataButton()
#         self.velScale()
#         self.canvas()

#     # print messages on the screen
#     def printMessage(self):
#         if (self.data):
#             print("Data is loaded and accessible from here (printMessage()).")
#         else:
#             print('No data loaded...')

#     ### OK button
#     def okButton(self):
#         self.okButton = tk.Button(self, text='Test', command=self.printMessage)
#         self.okButton.grid(column=1, row=1, sticky="nesw")

#     ### Quit button
#     def quitButton(self):
#         self.quitButton = tk.Button(self, text='Quit', command=self.confirmQuit)
#         self.quitButton.grid(column=1, row=2, sticky="nesw")
#     # confirm quitting
#     def confirmQuit(self):
#         answer = tkMessageBox.askyesno(title="App", message="Do you really want to quit?")
#         if (answer):
#             self.quit()

#     # Clear data from current session
#     def clearDataButton(self):
#         self.clearData = tk.Button(self, text='Clear data', command=self.confirmClearData)
#         self.clearData.grid(column=1, row=4, sticky="nesw")
#     # confirm clearing data
#     def confirmClearData(self):
#         answer = tkMessageBox.askyesno(title="App", message="Are you sure you want to clear the loaded data?")
#         if (answer):
#             self.data = None
#             tkMessageBox.showwarning(title="App", message="Data has been deleted.")

#     # Velocity scale
#     def velScale(self):
#         self.velVar = tk.StringVar()
#         velLabel = tk.Label(self, text="Scale value:", textvariable=self.velVar)
#         velLabel.grid(row=4, column=0, columnspan=2, sticky=tk.W+tk.E)
#         velScale = tk.Scale(self, from_=-500, to=+500, orient=tk.HORIZONTAL, resolution=20,
#                         sliderlength=20, showvalue=0,
#                         length=200, width=20,
#                         command=self.onVelScale)
#         velScale.grid(column=1, row=5, sticky="nesw")
        
#     # update velLabel
#     def onVelScale(self, val):
#         self.velVar.set("Scale value: {:+0.0f}".format(float(val)))

#     # Canvas
#     def canvas(self):
#         self.f = Figure(figsize=(4,2))
#         self.a = self.f.add_subplot(111)
#         self.a.plot([1,2,3,4,5,6,7,8],[5,6,1,3,8,9,3,5])

#         self.canvas = FigureCanvasTkAgg(self.f, master=self)
#         self.canvas.get_tk_widget().grid(column=2, row=1, rowspan=5, sticky="nesw")

# if __name__ == "__main__":
#     root = tk.Tk()
#     root.geometry("800x600+10+10")
#     root.resizable(0, 0)
#     MainApp(root).pack(side=tk.TOP)
#     root.mainloop()



# import tkinter as tk


# def on_configure(event):
#     # update scrollregion after starting 'mainloop'
#     # when all widgets are in canvas
#     canvas.configure(scrollregion=canvas.bbox('all'))


# root = tk.Tk()


#################################################################################### 

# # --- create canvas with scrollbar ---

# canvas = tk.Canvas(root)
# canvas.pack(side=tk.LEFT)

# scrollbar = tk.Scrollbar(root, command=canvas.yview)
# scrollbar.pack(side=tk.LEFT, fill='y')

# canvas.configure(yscrollcommand = scrollbar.set)

# # update scrollregion after starting 'mainloop'
# # when all widgets are in canvas
# canvas.bind('<Configure>', on_configure)

# # --- put frame in canvas ---

# frame = tk.Frame(canvas)
# canvas.create_window((0,0), window=frame, anchor='nw')

# # --- add widgets in frame ---

# l = tk.Label(frame, text="Hello", font="-size 50")
# l.pack()

# l = tk.Label(frame, text="World", font="-size 50")
# l.pack()

# l = tk.Label(frame, text="Test text 1\nTest text 2\nTest text 3\nTest text 4\nTest text 5\nTest text 6\nTest text 7\nTest text 8\nTest text 9", font="-size 20")
# l.pack()

# # --- start program ---

# root.mainloop()
