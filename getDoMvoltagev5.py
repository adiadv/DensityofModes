import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, sosfilt, welch, sosfiltfilt
from scipy.integrate import cumulative_trapezoid
from scipy.fftpack import fft, ifft

#UPDATE ME
#########################################################################################
rate = 500000
timespan = 2
nosamples = int(rate*timespan)

filenames = ['filename1.txt', 'filename2.txt', 'filename3.txt']
colors = ['r', 'g', 'b']
line_thickness = [1.5, 2.0, 2.5]
labels = ['File 1', 'File 2', 'File 3']
#########################################################################################

def timekeeper(n, deltat):
    return np.arange(n)*deltat

#########################################################################################

def getDoM(voltage_vec, rate_sc, duration_sc):

    #Voltage
    voltage_mean = np.mean(voltage_vec)
    voltage_vec -= voltage_mean

    #Velocity
    velocity_vec = cumulative_trapezoid(voltage_vec, initial=0)
    velocity_mean = np.mean(velocity_vec)
    velocity_vec -= velocity_mean

    #NYQ FILTER #2
    nyq = 0.5*rate_sc
    low = 47/nyq
    high = 68/nyq
    order = 2
    sos = butter(order, [low, high], btype= 'bandstop', output = 'sos')
    velocity_vec = sosfiltfilt(sos, velocity_vec)

    #padding
    velocity_vec_pad = np.zeros(len(velocity_vec))
    velocity_to_vac = np.concatenate((velocity_vec, velocity_vec_pad))

    #Velocity-Autocorrelation
    velocity_fft = fft(velocity_to_vac)
    velocity_fft *= np.conj(velocity_fft)
    velocity_fft_inv = ifft(velocity_fft)
    VAC_vec = velocity_fft_inv
    VAC_vec /= VAC_vec[0]

    #Density of Modes
    FFT_VAC_vec = fft(VAC_vec)
    DoM_vec = np.real(FFT_VAC_vec)
    P2 = np.abs(DoM_vec / int(rate*timespan))
    P1 = P2[:, int(rate*timespan) // 2 + 1]
    P1[1:-1] *= 2

    return voltage_vec, velocity_vec, VAC_vec, P1

#########################################################################################

#Initialize Plots
tv4 = timekeeper(nosamples, 1 / rate)
tv6 = timekeeper(2* nosamples, 1/rate)
fv1 = 0.5 * rate * np.arange(nosamples // 2 + 1) / nosamples

lenfv1 = len(fv1)

voltage_vec_outs = []
velocity_vec_outs = []
VAC_vec_outs = []
DoM_vec_outs = []
psdx_list = []
fvpsd_list = []

for fname in filenames:
    txtvoltage = np.loadtxt(fname)
    voltage_vec_in = txtvoltage[:,1]

    voltage_vec_out, velocity_vec_out, VAC_vec_out, DoM_vec_out = getDoM(voltage_vec_in, rate, timespan)

    voltage_vec_outs.append(voltage_vec_out)
    velocity_vec_outs.append(velocity_vec_out)
    VAC_vec_outs.append(VAC_vec_out)
    DoM_vec_outs.append(DoM_vec_out)

    xdft = fft(voltage_vec_out)
    xdft = xdft[:nosamples // 2 + 1]
    psdx = (1 / (rate * nosamples)) * np.abs(xdft) ** 2
    psdx[1:-1] *= 2
    fv3 = np.linspace(0, rate / 2, len(psdx))

    psdx_list.append(psdx)
    fvpsd_list.append(fv3)


#Voltage Plot
plt.figure()
for i, voltage_vec_out in enumerate(voltage_vec_outs):
    plt.plot(tv4, voltage_vec_out, color = colors[i], linewidth = line_thickness[i], label = labels[i])
plt.title('Voltage vs Time')
plt.ylabel('Voltage (V)')
plt.xlabel('Time (s)')
plt.show()

#Velocity Plot
plt.figure()
for i, velocity_vec_out in enumerate(velocity_vec_outs):
    plt.plot(tv4, velocity_vec_out, color = colors[i], linewidth = line_thickness[i], label = labels[i])
plt.title('Velocity v/s Time')
plt.ylabel('Velocity (au)')
plt.xlabel('Time (s)')
plt.show()

#VAC Plot
plt.figure()
for i, VAC_vec_out in enumerate(VAC_vec_outs):
    plt.plot(tv6, VAC_vec_out, color = colors[i], linewidth = line_thickness[i], label = labels[i])
plt.title('V.A.C. vs Delay Time')
plt.ylabel('Velocity Autocorrelation')
plt.xlabel('Time (s)')
plt.xlim([0, timespan / 2])
plt.ylim([-1, 1])
plt.show()

#DoM Plot
plt.figure()
for i, DoM_vec_out in enumerate(DoM_vec_outs):
    plt.loglog(fv1, DoM_vec_out, color = colors[i], linewidth = line_thickness[i], label = labels[i])
plt.title('Density of Modes')
plt.ylabel('Density of Modes')
plt.xlabel('Frequency (Hz)')
plt.xlim([0, rate / 2])
plt.show()

#Voltage Semilogy
xdft = fft(voltage_vec_out)
xdft = xdft[:nosamples//2 + 1]
psdx = (1/(rate * nosamples)) * np.abs(xdft)**2
psdx[1:-1] *= 2
fv3 = np.linspace(0,rate/2, len(psdx))

plt.figure()
for i, psdx in enumerate(psdx_list):
    plt.semilogy(fv3, psdx, color = colors[i], linewidth = line_thickness[i], label = labels[i])
plt.title('PSD vs f')
plt.ylabel('Power Spectral Density')
plt.xlabel('Frequency (Hz)')
plt.xlim([0, rate / 2])
plt.ylim([1e-16, 1])
plt.show()

#########################################################################################


