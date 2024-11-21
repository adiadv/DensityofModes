#########################################################################################
#Here we go again...
#My goal is to make this script cleaner, more efficient and less disastrous than the MATLAB version.

#This file will calculate
# 1. Voltages
# 2. Velocities
# 3. Velocity Autocorrelations
# 4. Density of Modes
#
# For a single device (The DoMinator) at a given rate and timespan (specified below).
# The file input must be a single column of voltages
#
# Here are the imports.
#########################################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, sosfilt, welch, sosfiltfilt
from scipy.integrate import cumulative_trapezoid
from scipy.fftpack import fft, ifft


#########################################################################################
#UPDATE THIS SECTION
rate = 100000000 #update this
timespan = 1 #update this
nosamples = int(rate*timespan)

filenames = ['rot5_dataonly.csv', 'rot4_dataonly.csv','rot3_dataonly.csv','rot2_dataonly.csv', 'rot1_dataonly.csv','rot0_dataonly.csv'] #add the files
colors = ['red','orange','gold','limegreen','b','m']#colors for the plots
line_thickness = [1.0,1.0,1.0,1.0,1.0,1.0] #leave this probably (imo)
labels = ['5% Strain', '4% Strain', '3% Strain', '2% Strain', '1% Strain', '0% Strain'] #labels for legends
#########################################################################################

def timekeeper(n, deltat): #as usual, leave the timekeeper be. If you decide to forgo the timekeeper, face the consequences on your own (or email me your choice)
    return np.arange(n)*deltat

#########################################################################################

def getDoM(voltage_vec, rate_sc, duration_sc): #NEW AND IMPROVED!

    #Voltage
    voltage_mean = np.mean(voltage_vec)
    voltage_vec -= voltage_mean #this never worked in matlab? dont ask questions.

    #Velocity
    velocity_vec = cumulative_trapezoid(voltage_vec, initial=0)
    velocity_mean = np.mean(velocity_vec)
    velocity_vec -= velocity_mean #subtract the mean again just to be sure because im not sure

    #NYQ FILTER #2- my choice (the other nuclear filters were just doing too much yk?)
    nyq = 0.5*rate_sc
    low = 55/nyq
    high = 65/nyq
    order = 2
    sos = butter(order, [low, high], btype= 'bandstop', output = 'sos')
    velocity_vec = sosfiltfilt(sos, velocity_vec)

    #padding
    velocity_vec_pad = np.zeros(len(velocity_vec))
    velocity_to_vac = np.concatenate((velocity_vec, velocity_vec_pad)) #slap on some zeros at the end

    #Velocity-Autocorrelation
    velocity_fft = fft(velocity_to_vac)
    velocity_fft *= np.conj(velocity_fft)
    velocity_fft_inv = ifft(velocity_fft) #oh yeah
    VAC_vec = np.real(velocity_fft_inv) #this was the prooobllemeeemmem i seeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
    VAC_vec /= VAC_vec[0] #normalize


    #Density of Modes
    '''
    FFT_VAC_vec = fft(VAC_vec)
    DoM_vec = np.real(FFT_VAC_vec)
    P2 = np.abs(DoM_vec / int(rate*timespan))
    P1 = P2[:int(rate*timespan) // 2 + 1]
    P1[1:-1] *= 2
    '''

    fv_plot, DoM_vec = welch(VAC_vec, fs = rate_sc, nperseg = len(VAC_vec)//8, scaling = 'density')
    DoM_vec = np.real(DoM_vec)
    #oops where is the factor of 2 (me and fourier space)

    return voltage_vec, velocity_vec, VAC_vec, fv_plot, DoM_vec
    #I think I ended up finding the factor of two but not sure. If there is a problem just divide everything you need by 2 until it works

#########################################################################################

#Initialize Plots
tv4 = timekeeper(nosamples, 1 / rate) #I dont think we need tv1, tv2, tv5, tv7 so ill keep what I need.
tv6 = timekeeper(2* nosamples, 1/rate)
fv1 = 0.5 * rate * np.arange(nosamples // 2 + 1) / nosamples

lenfv1 = len(fv1)

voltage_vec_outs = [] #INITIALIZE
velocity_vec_outs = [] #ARRAYS
VAC_vec_outs = [] #FOR
fv_plot_outs = [] #FILLING
DoM_vec_outs = [] #WITH
psdx_list = [] #GETDOM
fvpsd_list = [] #DATA :)

for fname in filenames:
    txtvoltage = pd.read_csv(fname, header = None) #you can also use np.loadtxt
    voltage_vec_in = txtvoltage.iloc[:,1].to_numpy() #lock in and edit if needed based on your file type (dont do that actually)

    voltage_vec_out, velocity_vec_out, VAC_vec_out, fv_plot_out, DoM_vec_out = getDoM(voltage_vec_in, rate, timespan)

    voltage_vec_outs.append(voltage_vec_out) #APPEND
    velocity_vec_outs.append(velocity_vec_out) #GETDOM
    VAC_vec_outs.append(VAC_vec_out) #OUTS
    fv_plot_outs.append(fv_plot_out) #FOR
    DoM_vec_outs.append(DoM_vec_out) #PLOTTING

    xdft = fft(voltage_vec_out) #honestly i cant rememebr what this is but Im taking it from the matlab script (its the voltage semilogy stuffs)
    xdft = xdft[:nosamples // 2 + 1]
    psdx = (1 / (rate * nosamples)) * np.abs(xdft) ** 2
    psdx[1:-1] *= 2
    fv3 = np.linspace(0, rate / 2, len(psdx))

    psdx_list.append(psdx) #:)
    fvpsd_list.append(fv3) #:(


#Voltage Plot
plt.figure()
for i, voltage_vec_out in enumerate(voltage_vec_outs):
    plt.plot(tv4, voltage_vec_out, color = colors[i], linewidth = line_thickness[i], label = labels[i]) #now you see a random spike and realize you have to go retake data because the DoMinator fell asleep :)
plt.title('Voltage vs Time')
plt.ylabel('Voltage (V)')
plt.xlabel('Time (s)')
plt.legend(loc = 'best')
plt.show()

#Velocity Plot
plt.figure()
for i, velocity_vec_out in enumerate(velocity_vec_outs):
    plt.plot(tv4, velocity_vec_out, color = colors[i], linewidth = line_thickness[i], label = labels[i]) #now you see a giant velocity value and realize you need to go retake data
plt.title('Velocity v/s Time')
plt.ylabel('Velocity (au)')
plt.xlabel('Time (s)')
plt.legend(loc = 'best')
plt.show()

#VAC Plot
plt.figure()
for i, VAC_vec_out in enumerate(VAC_vec_outs):
    plt.plot(tv6, VAC_vec_out, color = colors[i], linewidth = line_thickness[i], label = labels[i]) #Then you see a flat line and realize that something else went wrong (wtf?) so you need to retake data
plt.title('V.A.C. vs Delay Time')
plt.ylabel('Velocity Autocorrelation')
plt.xlabel('Time (s)')
plt.xlim([0, timespan / 2])
plt.ylim([-1, 1])
plt.legend(loc = 'best')
plt.show()

#DoM Plot
plt.figure()
for i in range(len(DoM_vec_outs)):
    plt.loglog(fv_plot_outs[i], 10**(10*i) * DoM_vec_outs[i], color = colors[i], linewidth = line_thickness[i], label = labels[i]) #now the DoM is plotted with some 60hz noise but you edit it and fix that and then realize that its sloped downwards? Why? is my system cooked? #Debye scaling who?
plt.title('Density of Modes')
plt.ylabel('Density of Modes (au)')
plt.xlabel('Frequency (Hz)')
plt.xlim([1, rate / 2])
#plt.xlim([8000, 800000])
plt.legend(loc = 'best')
plt.show()


#Voltage Semilogy
fv3 = np.linspace(0,rate/2, len(psdx))

plt.figure()
for i in range(len(psdx_list)):
    plt.loglog(fvpsd_list[i], psdx_list[i], color = colors[i], linewidth = line_thickness[i], label = labels[i]) #lets be honest you dont check this if you dont have to because its after the DoM plot and thats all you care about isnt it
plt.title('PSD vs f')
plt.ylabel('Power Spectral Density')
plt.xlabel('Frequency (Hz)')
plt.xlim([0, rate / 2])
plt.ylim([1e-16, 1])
plt.legend(loc = 'best')
plt.show()


#########################################################################################


