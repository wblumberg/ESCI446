#get_ipython().getoutput("/usr/bin/env python")
# File: FFT_demonstrator.py

'''This Tkinter application demonstrates the recovery of a digitized function by summing its spectral
component.'''

# Written by Alex DeCaria
# Latest version date:  3/1/2019

# NOTE:  Make sure PyLab is disabled

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
class FFT_demonstrator:

    def __init__(self, master):

        self.n = 101   #  Number of data points  (KEEP ODDget_ipython().getoutput(")")
        self.factor = 10 #  Used for scaling plots
        self.x = np.arange(0,self.factor*(self.n-1)+1)/float(self.factor)
        self.signal = np.zeros_like(self.x)
        self.pre_y = np.zeros_like(self.signal)
        
        self.main_figure = fig.Figure(figsize=(4,1.2))  #  this creates the main figure for the plot
        self.real_figure = fig.Figure(figsize=(4,1.2))  #  this creates the main figure for the plot
        self.imaginary_figure = fig.Figure(figsize=(4,1.2))  #  this creates the main figure for the plot
        self.pre_figure = fig.Figure(figsize=(2,1)) # this creates preview figure for the plot

        master.title('FFT Demonstrator')
        master.resizable(0,0)

        # Creates an outer frame within the main window
        self.outer_frame = Frame(master)
        self.outer_frame.pack(side = LEFT, expand=YES, fill=BOTH)

        # Creates a frame for figures
        self.figure_frame = Frame(self.outer_frame)
        self.figure_frame.pack(side = LEFT, expand=YES, fill=X)

        # Creates a label frame for main figure
        self.mainfig_frame = LabelFrame(self.figure_frame, text = 'Signal')
        self.mainfig_frame.pack(side = TOP, expand=YES, fill=X)

        # Creates a canvas for displaying main figure
        self.main_canvas = tkagg.FigureCanvasTkAgg(self.main_figure,self.mainfig_frame)
        self.main_canvas.draw()
        self.main_canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        # Creates a label frame for real spectral figure
        self.real_frame = LabelFrame(self.figure_frame, text = 'Real Spectral Coefficients')
        self.real_frame.pack(side = TOP, expand=YES, fill=X)

        # Creates a canvas for displaying real spectral figure
        self.real_canvas = tkagg.FigureCanvasTkAgg(self.real_figure,self.real_frame)
        self.real_canvas.draw()
        self.real_canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        # Creates a label frame for imaginary spectral figure
        self.imaginary_frame = LabelFrame(self.figure_frame, text = 'Imaginary Spectral Coefficients')
        self.imaginary_frame.pack(side = TOP, expand=YES, fill=X)

        # Creates a canvas for displaying imaginary spectral figure
        self.imaginary_canvas = tkagg.FigureCanvasTkAgg(self.imaginary_figure,self.imaginary_frame)
        self.imaginary_canvas.draw()
        self.imaginary_canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        # Creates a frame for holding preview and controls
        self.control_frame = LabelFrame(self.outer_frame, text = 'Wave Parameters')
        self.control_frame.pack(side = RIGHT, expand=YES, fill=X)

        # Creates a canvas for displaying preview figure
        self.pre_canvas = tkagg.FigureCanvasTkAgg(self.pre_figure,self.control_frame)
        self.pre_canvas.draw()
        self.pre_canvas.get_tk_widget().pack(side=TOP, fill=X, expand=YES)

        # Create frame for holding sliders
        self.slider_frame = Frame(self.control_frame)
        self.slider_frame.pack(side = TOP)

        # Create frame and widgets for choosing signal wavelength
        self.length_frame = LabelFrame(self.slider_frame, text = 'Cycles')
        self.length_frame.pack(side = LEFT, padx = 5, pady = 5)
        self.length_scale = Scale(self.length_frame, orient = HORIZONTAL, from_=1, to=50,
                        bigincrement = 10, command = self.create_preview)
        self.length_scale.set(2)
        self.length_scale.pack(side = LEFT)

        # Create frame and widgets for choosing signal amplitude
        self.amplitude_frame = LabelFrame(self.slider_frame, text = 'Amplitude')
        self.amplitude_frame.pack(side = LEFT, padx = 5, pady = 5)
        self.amplitude_scale = Scale(self.amplitude_frame, orient = HORIZONTAL, from_=1, to=10,
                        bigincrement = 1, command = self.create_preview)
        self.amplitude_scale.set(5)
        self.amplitude_scale.pack(side = LEFT)

        # Create frame and widgets for choosing signal phase
        self.phase_frame = LabelFrame(self.slider_frame, text = 'Phase (degress)')
        self.phase_frame.pack(side = LEFT, padx = 5, pady = 5)
        self.phase_scale = Scale(self.phase_frame, orient = HORIZONTAL, from_=-180, to=180,
                        bigincrement = 30, command = self.create_preview)
        self.phase_scale.set(0)
        self.phase_scale.pack(side = LEFT)
     
        # Create frame for holding buttons
        self.button_frame = Frame(self.control_frame)
        self.button_frame.pack(side = TOP)

        # Create button for adding wave to signal and plotting
        self.plot_button = Button(self.button_frame, text = 'Add Wave',
                                     command = self.add_to_signal)
        self.plot_button.pack(side = LEFT)

        # Create button for resetting plot
        self.reset_button = Button(self.button_frame, text = 'RESET',
                                     command = self.reset)
        self.reset_button.pack(side = LEFT)

     
    def create_preview(self, val=0):
        self.pre_y = np.ndarray(self.factor*(self.n-1)+1, dtype = np.float64)
        amp = float(self.amplitude_scale.get())           # real amplitude
        p = np.radians(float(self.phase_scale.get()))   #  phase in radians
        l = (self.n-1)/float(self.length_scale.get())
        k = 2*np.pi/float(l)      # angular wavenumber
        for i, j in enumerate(self.x):
            arg = np.complex(0, k*j - p)
            self.pre_y[i] = np.real(amp*np.exp(arg))
        self.pre_figure.clf()
        a = self.pre_figure.add_axes([0.1, 0.25, 0.8, 0.6])
        a.plot(self.x,self.pre_y)
        a.set_ylim(-10, 10)
        a.set_xlim(self.x[0],self.x[-1])
        self.pre_canvas.draw()
        
    def add_to_signal(self, val=0):
        self.signal += self.pre_y
        self.main_figure.clf()
        a = self.main_figure.add_axes([0.1, 0.2, 0.8, 0.65])
        a.plot(self.x, self.signal)
        a.set_xlim(self.x[0],self.x[-1])
        self.main_canvas.draw()
        self.spectral_update()

    def spectral_update(self):
        s = np.zeros(self.n-1, dtype = np.float64)  # coarse signal
        f = np.fft.fftfreq(len(s))  #  frequencies
        f = np.fft.fftshift(f) #  shifts zero freq to center
        f = np.append(f,-f[0])  # adds positive Nyquist f to end of freqs.
        
        for i in range(0,self.n-1):
            s[i] = self.signal[i*self.factor]
        Sn = np.fft.fft(s)/float(len(s))  # fourier coefficients
        Sn = np.fft.fftshift(Sn)  #  shifts zero freq to center
        Sn = np.append(Sn, np.conj(Sn[0]))  # adds complex conjugate of first element to end of array
        maxreal = np.max(np.abs(np.real(Sn)))
        maximag = np.max(np.abs(np.imag(Sn)))
        maxboth = max(maxreal, maximag)
        if maxboth == 0:
            maxboth = 1.0
        
        self.real_figure.clf()
        a = self.real_figure.add_axes([0.1,0.2,0.8,0.7])
        a.bar(f, np.real(Sn),color = 'k', width = 0.005, lw = 0, align = 'center')
        mn,mx = a.set_xlim()
        a.hlines(0, mn, mx, color='brown', linestyle = ':') 
        a.set_ylim(-maxboth, maxboth)
        self.real_canvas.draw()

        self.imaginary_figure.clf()
        b = self.imaginary_figure.add_axes([0.1,0.2,0.8,0.7])
        b.bar(f, np.imag(Sn),color = 'k', width = 0.005, lw = 0, align = 'center')
        mn,mx = b.set_xlim()
        b.hlines(0, mn, mx, color='brown', linestyle = ':') 
        b.set_ylim(-maxboth, maxboth)
        self.imaginary_canvas.draw()

    def reset(self, val=0):
        self.signal = 0*self.signal
        self.main_figure.clf()
        a = self.main_figure.add_subplot(111)
        a.plot(self.x, self.signal)
        self.main_canvas.draw()
        self.spectral_update()
        
        
    
# ------------ END OF APP DEFINITION -------------------------------------------------------

# Calls Tk library, instantiates application, and enters the main loop
root = Tk()
app = FFT_demonstrator(root)
root.mainloop()








