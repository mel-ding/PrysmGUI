# Starter Code from: https://pythonprogramming.net/how-to-embed-matplotlib-graph-tkinter-gui/

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import tkinter as tk
from tkinter import ttk

LARGE_FONT= ("Verdana", 25)


class Grapher(tk.Tk):
    # Frame handling logic here: 
    # https://stackoverflow.com/questions/7546050/switch-between-two-frames-in-tkinter

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        # the container is where we'll stack a bunch of frames on top of each 
        # other, then the one we want visible will be raised above the others
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (StartPage, PageOne, PageTwo, PageThree):
            # put all of the pages in the same location;
            # the one on the top of the stacking order
            # will be the one that is visible.
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()


class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        label = tk.Label(self, text="Start Page", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button = ttk.Button(self, text="Visit Page 1",
                            command=lambda: controller.show_frame(PageOne))
        button.pack()

        button2 = ttk.Button(self, text="Visit Page 2",
                            command=lambda: controller.show_frame(PageTwo))
        button2.pack()

        button3 = ttk.Button(self, text="Graph Page",
                            command=lambda: controller.show_frame(PageThree))
        button3.pack()

        # widget.bind(event, handler)
        # read more here: http://effbot.org/tkinterbook/tkinter-events-and-bindings.htm
        # events are strings: <modifier>
        self.bind("<Button-1>", self.left_click)

    def left_click(self, event):
        tk.Label(self, text="Left Click!").pack()


class PageOne(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Page One!!!", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()

        button2 = ttk.Button(self, text="Page Two",
                            command=lambda: controller.show_frame(PageTwo))
        button2.pack()

        button3 = ttk.Button(self, text="Event Clicker",
                            command=lambda: controller.show_frame(StartPage))
        button3.pack()



class PageTwo(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Page Two!!!", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()

        button2 = ttk.Button(self, text="Page One",
                            command=lambda: controller.show_frame(PageOne))
        button2.pack()


class PageThree(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Graph Page!", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self, text="Back to Home",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()

        # Coefficient Entry Boxes
        tk.Label(self, text = "x^2 coefficient").pack()
        self.aEntry = tk.Entry(self)
        self.aEntry.pack()

        tk.Label(self, text = "x coefficient").pack()
        self.bEntry = tk.Entry(self)
        self.bEntry.pack()
        
        tk.Label(self, text = "constant coefficient").pack()
        self.cEntry = tk.Entry(self)
        self.cEntry.pack()

        # Creates graph
        graphButton = tk.Button(self, text="Graph", command=self.quadraticVals)
        graphButton.pack()

        # Clears graph
        clearButton = tk.Button(self, text="Clear", command=self.clear)
        clearButton.pack()

        f = Figure(figsize=(5,5), dpi=100)
        self.a = f.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(f, self)
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

    """
    Takes the given entry constants and creates a quadractic graph. 
    """
    def quadraticVals(self): 
        # clear whatever is currently on the canvas 
        self.a.clear()
        # get the entry fields 
        aCoeff = int(self.aEntry.get())
        bCoeff = int(self.bEntry.get())
        cCoeff = int(self.cEntry.get())
        # create 1000 equally spaced points between -10 and 10 
        xVals = np.linspace(-20, 20, 1000)
        # calculate y value for each element in x vector 
        yVals = (aCoeff * xVals**2) + (bCoeff * xVals) + cCoeff
        # plot the graph 
        self.a.plot(xVals, yVals)
        self.canvas.draw()

    """
    Clears the content in the graph. 
    """
    def clear(self):
        self.a.clear()
        self.canvas.draw()

app = Grapher()
app.mainloop()