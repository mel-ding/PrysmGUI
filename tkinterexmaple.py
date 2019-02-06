# from Tkinter import *
#
#
# window  = Tk()
# window.title("Visualize Quadratic Functions")
# lbl = Label(window, text = "Hello", font= ("Arial Bold", 50))
# lbl.grid(column = 0, row=0)
# window.geometry('500x500')
# def clicked():
#     lbl.configure(text= "button was clicked!")
# txt = Entry(window, width=10)
# mainframe = ttk.Frame(window, padding="3 3 12 12")
# mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
# window.columnconfigure(0, weight=1)
# window.rowconfigure(0, weight=1)
# btn = Button(window, text="Click Me", command = clicked)
# btn.grid(column=1, row=0)
# window.mainloop()
##########################################################################################################
# from tkinter import *
# from tkinter import ttk
# def calculate(*args):
#     try:
#         value1 = float(coff_quadratic.get())
#         value2 = float(coff_linear.get())
#         value3 = float(constant.get())
#         meters.set((0.3048 * value1 * 10000.0 + 0.5)/10000.0)
#     except ValueError:
#         pass
#
# root = Tk()
# root.title("Quadratic visualizer")
#
# mainframe = ttk.Frame(root, padding="3 3 12 12")
# mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
# root.columnconfigure(0, weight=1)
# root.rowconfigure(0, weight=1)
#
# coff_quadratic = StringVar()
# coff_linear = StringVar()
# constant = StringVar()
# meters = StringVar()
#
# quadratic_entry = ttk.Entry(mainframe, width=7, textvariable=coff_quadratic)
# linear_entry = ttk.Entry(mainframe, width=7, textvariable=coff_linear)
# constant_entry = ttk.Entry(mainframe, width=7, textvariable=constant)
#
# quadratic_entry.grid(column=2, row=1, sticky=(W, E))
# linear_entry.grid(column=2, row=1, sticky=(W, E))
# constant_entry.grid(column=2, row=1, sticky=(W, E))
#
# ttk.Label(mainframe, textvariable=meters).grid(column=2, row=2, sticky=(W, E))
# ttk.Button(mainframe, text="Calculate", command=calculate).grid(column=3, row=3, sticky=W)
#
# ttk.Label(mainframe, text="feet").grid(column=3, row=1, sticky=W)
# ttk.Label(mainframe, text="is equivalent to").grid(column=1, row=2, sticky=E)
# ttk.Label(mainframe, text="meters").grid(column=3, row=2, sticky=W)
#
# for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)
#
# quadratic_entry.focus()
# root.bind('<Return>', calculate)
#
# root.mainloop()
################################################################################################
import matplotlib as mpl
import numpy as np
import sys
if sys.version_info[0] < 3:
    import Tkinter as tk
else:
    import tkinter as tk
import matplotlib.backends.tkagg as tkagg
from matplotlib.backends.backend_agg import FigureCanvasAgg


def draw_figure(canvas, figure, loc=(0, 0)):
    """ Draw a matplotlib figure onto a Tk canvas

    loc: location of top-left corner of figure on canvas in pixels.
    Inspired by matplotlib source: lib/matplotlib/backends/backend_tkagg.py
    """
    figure_canvas_agg = FigureCanvasAgg(figure)
    figure_canvas_agg.draw()
    figure_x, figure_y, figure_w, figure_h = figure.bbox.bounds
    figure_w, figure_h = int(figure_w), int(figure_h)
    photo = tk.PhotoImage(master=canvas, width=figure_w, height=figure_h)

    # Position: convert from top-left anchor to center anchor
    canvas.create_image(loc[0] + figure_w/2, loc[1] + figure_h/2, image=photo)

    # Unfortunately, there's no accessor for the pointer to the native renderer
    tkagg.blit(photo, figure_canvas_agg.get_renderer()._renderer, colormode=2)

    # Return a handle which contains a reference to the photo object
    # which must be kept live or else the picture disappears
    return photo

# Create a canvas
w, h = 300, 200
window = tk.Tk()
window.title("A figure in a canvas")
canvas = tk.Canvas(window, width=w, height=h)
canvas.pack()

# Generate some example data
X = np.linspace(0, 2 * np.pi, 50)
Y = np.sin(X)

# Create the figure we desire to add to an existing canvas
fig = mpl.figure.Figure(figsize=(2, 1))
ax = fig.add_axes([0, 0, 1, 1])
ax.plot(X, Y)

# Keep this handle alive, or else figure will disappear
fig_x, fig_y = 100, 100
fig_photo = draw_figure(canvas, fig, loc=(fig_x, fig_y))
fig_w, fig_h = fig_photo.width(), fig_photo.height()

# Add more elements to the canvas, potentially on top of the figure
canvas.create_line(200, 50, fig_x + fig_w / 2, fig_y + fig_h / 2)
canvas.create_text(200, 50, text="Zero-crossing", anchor="s")

# Let Tk take over
tk.mainloop()
