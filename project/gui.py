import tkinter as tk
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ternary
import time
import PIL as pl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def get_input():
    # Get the values from the entry widgets and store them in variables
    value1 = float(entry1.get())
    value2 = float(entry2.get())
    value3 = float(entry3.get())
    value4 = float(entry4.get())
    value5 = float(entry5.get())
    value6 = float(entry6.get())
    toggle_value = float(entry7.get())
    # Get the value of the toggle button and store it in a variable
    output_label.config(text="Wait for 40 seconds to get the output printed")

    if toggle_value == 1:
        os.system('cmd /c "matlab -nosplash -nodesktop -r "func_counter({},{},{},{},{});exit;" "'.format(value1, value2, value3, value4, value6))
        time.sleep(40)
        data = pd.read_excel('out_counter.xlsx', header=None)
        data = data.to_numpy()
        percent_removal = data[0,0]
        stage = int(data[0,1])
        plot = plt.imread('plot_counter.png')
        fig,ax = plt.subplots(1)
        ax.imshow(plot)
        plt.axis('off')
        flow = lambda x: 'Cross Current' if toggle_value==0 else 'Counter Current'
        output_label.config(text="Percentage Removal - {}\nFlow Type : {}".format(percent_removal, flow(toggle_value)))
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack()

    else:
        os.system('cmd /c "matlab -nosplash -nodesktop -r "func_cross({},{},{},{},{});exit;" "'.format(value1, value2, value3, value4, value5))
        time.sleep(40)
        data = pd.read_excel('out_cross.xlsx', header=None)
        data = data.to_numpy()
        percent_removal = data[0,0]
        stage = int(data[0,1])
        plot = plt.imread('plot_cross.png')
        fig,ax = plt.subplots(1)
        ax.imshow(plot)
        plt.axis('off')
        flow = lambda x: 'Cross Current' if toggle_value==0 else 'Counter Current'
        output_label.config(text="Percentage Removal - {}\nFlow Type : {}".format(percent_removal, flow(toggle_value)))
        canvas = FigureCanvasTkAgg(fig, master=root)
        canvas.draw()
        canvas.get_tk_widget().pack()


# Create the main window
root = tk.Tk()
root.geometry("720x1680")
# Set the window title
root.title("Mass Transfer Project Group 2")

# Create label widgets for each input
label1 = tk.Label(root, text="F :")
label1.pack()
entry1 = tk.Entry(root)
entry1.pack()

label2 = tk.Label(root, text="S :")
label2.pack()
entry2 = tk.Entry(root)
entry2.pack()

label3 = tk.Label(root, text="xcf :")
label3.pack()
entry3 = tk.Entry(root)
entry3.pack()

label4 = tk.Label(root, text="xcs :")
label4.pack()
entry4 = tk.Entry(root)
entry4.pack()

label5 = tk.Label(root, text="stage :(0 if counter current)")
label5.pack()
entry5 = tk.Entry(root)
entry5.pack()

label6 = tk.Label(root, text="xcrn :(0 if cross current)")
label6.pack()
entry6 = tk.Entry(root)
entry6.pack()

# Create a toggle button widget
label7 = tk.Label(root, text="For Cross : 0 & Counter : 1 ")
label7.pack()
entry7 = tk.Entry(root)
entry7.pack()

# Create a button widget
button = tk.Button(root, text="Submit", command=get_input)
button.pack()

# Create a label widget for the output
output_label = tk.Label(root, text="")
output_label.pack()

# Start the main event loop
root.mainloop()