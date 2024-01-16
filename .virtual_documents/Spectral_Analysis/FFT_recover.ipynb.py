#get_ipython().getoutput("/usr/bin/env python")
# File: FFT_recover.py

'''This Tkinter application demonstrates the recovery of a digitized function by summing its spectral
component.'''

# Written by Alex DeCaria
# Latest version date:  3/26/2019

# NOTE!  Make sure PyLab is disabledget_ipython().getoutput("!")

import sys
if sys.version_info[0] get_ipython().getoutput("= 3:")
    from Tkinter import *    # if Python 2.X
else:
    from tkinter import *    # if Python 3.X
    
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure as fig
import matplotlib.backends.backend_tkagg as tkagg

# Define class for GUI Application
class FFT_recover:

    def __init__(self, master):

        self.n = 101   #  Number of data points
        self.x = np.arange(0,self.n) # array for x values
        self.signal = np.zeros_like(self.x, dtype = np.float64)
        self.S = np.zeros_like(self.x, dtype = np.complex)
        
        self.f = fig.Figure(figsize=(4,3))  #  this creates the figure for the plot

        master.title('FFT Recovery')
        #master.resizable(0,0)

        # Creates an outer frame within the main window
        self.outer_frame = Frame(master)
        self.outer_frame.pack(expand=YES, fill=BOTH)

        # Creates a frame for displaying plots
        self.canvas = tkagg.FigureCanvasTkAgg(self.f,self.outer_frame)
        #self.canva.show()  # This was deprecated
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=YES)

        # Create container frame for radio buttons and slider widgets
        self.widget_frame = Frame(self.outer_frame)
        self.widget_frame.pack(side = BOTTOM, expand=YES)

        # Create frame and widgets for choosing signal shape
        self.shape_frame = LabelFrame(self.widget_frame, text = 'Signal Shape')
        self.shape_frame.pack(side = LEFT, fill = BOTH, padx = 5, pady = 5, expand=YES)
        self.shape = IntVar()  #  Variable to hold state.  0 = square, 1 = bell
        self.shape.set(0)
        #print(self.shape.get())
        Radiobutton(self.shape_frame, text = 'Square', variable = self.shape,
                    value = 0, command = self.calculations).pack(side = LEFT)
        Radiobutton(self.shape_frame, text = 'Bell', variable = self.shape,
                    value = 1, command = self.calculations).pack(side = LEFT, expand=YES)

        # Create frame and widgets for choosing signal width
        self.width_frame = LabelFrame(self.widget_frame, text = 'Signal Width')
        self.width_frame.pack(side = LEFT, fill = BOTH, padx = 5, pady = 5, expand=YES)
        self.width_scale = Scale(self.width_frame, orient = HORIZONTAL, from_=1, to=90,
                        bigincrement = 10, command = self.calculations)
        self.width_scale.set(20)
        self.width_scale.pack(fill = BOTH, expand=YES)

        # Create frame and widgets for choosing components
        self.component_frame = LabelFrame(self.widget_frame, text = 'Components')
        self.component_frame.pack(side = LEFT, fill = BOTH, padx = 5, pady = 5, expand=YES)
        self.component_scale = Scale(self.component_frame, orient = HORIZONTAL,
                                     from_ = 0, to = int(self.n/2), command = self.plot_figure)
        self.component_scale.set(0)
        self.component_scale.pack(fill = BOTH, expand=YES)
        
    def calculations(self, val = 0):
        # Need value because some widgets automatically send state
        w = self.width_scale.get()  #  Width of signal
        # Generate signal
        if self.shape.get() == 0:
            self.signal = np.zeros_like(self.x, dtype = np.float64)
            self.signal[int((self.n-w)/2):int((self.n+w)/2)] = 1.0
        else:
            for i in range(0,len(self.x)):
                self.signal[i] = np.exp(-((float(i)-self.n/2)/(w/5.0))**2)
        self.component_scale.set(0)  # Reset # of FFT components to plot to zero
        # Calculate FFT
        self.S = np.fft.fft(self.signal)/self.n  #  FFT coefficients
        self.k = np.fft.fftfreq(self.n)  #  FFT frequencies
        self.plot_figure()            

    def plot_figure(self, val = 0):
        # clear and recreate figure
        self.f.clf() 
        self.a = self.f.add_subplot(111)
        self.a.plot(self.x,self.signal, ':')  # plots signal
        self.a.set_ylim(-0.25, 1.25)
        self.a.set_xlim(self.x[0],self.x[-1])

        # Generate recovered signal with specified number of components
        n = self.component_scale.get() #  number of fft components to plot
        recovered = np.ones_like(self.signal)*np.real(self.S[0])  # mean of signal
        for i in range(1,n+1):
            recovered += np.real(self.S[i])*np.cos(2*np.pi*self.k[i]*self.x) + \
                         np.real(self.S[-i])*np.cos(2*np.pi*self.k[-i]*self.x) - \
                         np.imag(self.S[i])*np.sin(2*np.pi*self.k[i]*self.x) - \
                         np.imag(self.S[-i])*np.sin(2*np.pi*self.k[-i]*self.x)
                                 
        self.a.plot(self.x, recovered, '-')

        self.canvas.draw()
        
# ------------ END OF APP DEFINITION -------------------------------------------------------

# Calls Tk library, instantiates application, and enters the main loop
root = Tk()
app = FFT_recover(root)
root.mainloop()




