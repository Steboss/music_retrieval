import scipy.io.wavfile
import matplotlib.pyplot  as plt
from mpl_toolkits.mplot3d import Axes3D
import pydub
import stft
import numpy  as np
import math
import sys,os
import seaborn as sbn
sbn.set_style("whitegrid")
plt.rcParams["xtick.labelsize"]=20
plt.rcParams["ytick.labelsize"]=20


def readmp3(f, normalized=False):
    """MP3 to numpy array  16bit only"""
    a = pydub.AudioSegment.from_mp3(f)
    y = np.array(a.get_array_of_samples())
    if a.channels == 2:
        y = y.reshape((-1, 2))
    if normalized:
        return a.frame_rate, np.float32(y) / 2**15
    else:
        return a.frame_rate, y

def writemp3(f, sr, x, normalized=False):
    """numpy array to MP3"""
    channels = 2 if (x.ndim == 2 and x.shape[1] == 2) else 1
    if normalized:  # normalized array - each item should be a float in [-1, 1)
        y = np.int16(x * 2 ** 15)
    else:
        y = np.int16(x)
    song = pydub.AudioSegment(y.tobytes(), frame_rate=sr, sample_width=2, channels=channels)
    song.export(f, format="mp3", bitrate="320k")


def padding(x, windowSize):
    #iinitialy add the boundaries of the signal,
    #so create an array of windowSize/2 zeroes that will be added ad the top and
    #bottom of the signal
    #final singal [ 0,0,0,0,0,.....  x ..... 0,0,0,0,0,0] as many 0 as windowSize/2
    #on both sides

    zeros = np.zeros(int(windowSize/2), dtype=x.dtype)
    print(x.shape)
    print(zeros.shape)
    x = np.concatenate((zeros, x, zeros), axis=-1)
    #then
    # Pad to integer number of windowed segments
    # I.e make x.shape[-1] = nperseg + (nseg-1)*nstep, with integer nseg
    nadd = (-(x.shape[-1]-windowSize) % int(windowSize/2)) % int(windowSize)
    zeros_shape = list(x.shape[:-1]) + [int(nadd)]
    x = np.concatenate((x, np.zeros(zeros_shape)), axis=-1)
    return x


#MAIN#
#read the input file
#define the extension
inputFile = "input_songs/equation.mp3"
extension = inputFile.split(".")[-1]
if extension =="wav":
    rate, audData = scipy.io.wavfile.read(inputFile)
else:
    rate, audData = readmp3(inputFile)
#extract the info
length = int(len(audData))#/265152#int(len(audData))#audData.shape[0] # this is the number of samples,
#if you divide length by the rateyou get the length in seconds /rate
#print(length, length/rate)
#sys.exit(-1)
channel1 = audData[:,0][int(length/2):(length)]
print(channel1.shape)
print(length)
#scipy.io.wavfile.write("study.wav",rate,channel1)
print("Total number of elements in the signal %d"%length);
#convert in double format
channel1 = np.double(channel1)
#channel2 = audData[:,1]
windowSize = 4096 #length of the window to analyse with STFT
hopSize = 4096  #hopsize between windows

#pad the signal based on the windowSize
print("Padding signal...")
signal = padding(channel1, windowSize)
print("Dimension of the signal after boundary and padding %d\n"% len(signal))

#UNCOMMENT if you want to analyse just a piece of the entire wav
#save_wav = audData[:,0][0:windowSize]#[length:(length*2)]
#scipy.io.wavfile.write("study.wav",rate,save_wav)
#window_start = 3*12288
#window_end = window_start + 12288

#compute the STFT
print("STFT...")
#signal = [1,2,3,4,5,6,7,8,9,10]
#signal = np.asarray(signal)
#signal = padding(signal,8)
#print(signal)
#rate = 1
#windowSize = 8
#hopSize = 8

#create the frequencies
freq_spacing = 1.0/(rate) #the spacing
val_freq = 1.0/(windowSize*freq_spacing)#the total number of freq
N_val = math.floor(windowSize/2) + 1#
print("Total number of freqs %d\n", N_val)
frequencies = []
for i in range(0, N_val):
    frequencies.append(i*val_freq)

magnitude = stft.play(signal, rate , windowSize, hopSize)
print(max(frequencies))
#print(magnitude.shape)
print("Roll magnitude values")
magnitude = np.rollaxis(magnitude, -1, -2)
#print(magnitude.shape)
start_time = int(windowSize/2)
stop_time  = int(magnitude.shape[-1] - windowSize/2 +1)

time = np.arange(windowSize/2, signal.shape[-1] - windowSize/2 + 1,
                 windowSize - (windowSize/2))/float(rate)

time -= (windowSize/2)/ rate
#for i in range(0,10):#
#    print(magnitude[i])
#PLOT
#x-axis for  the data
'''
x_axis_ref = np.linspace(0, len(magnitude),len(magnitude))
#colors andfigure
colors= sbn.color_palette()
col_ref = sbn.color_palette("cubehelix", 8)
fig = plt.figure(figsize=[15,10])
ax = fig.add_subplot(111)
#axis plot
ax.plot(x_axis_ref,np.abs(magnitude))
ax.xaxis.set_tick_params(labelsize=30)
ax.yaxis.set_tick_params(labelsize=30)
ax.set_ylabel(r"Freq",fontsize=30)
ax.set_xlabel(r"ticks",fontsize=30)
#lgd = ax.legend(loc="best", fontsize=30)#ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          #fancybox=True, shadow=True, ncol=3,fontsize=30)

plt.tight_layout()
plt.savefig("power_spectrum.pdf")#,bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.clf()
'''
#plot frequenceis
print("plotting frequencies...")
fig = plt.figure(figsize=[15,10])
ax = fig.add_subplot(111)
#ax.imshow(np.abs(magnitude))
print("Computing magnitude in dB...")
magnitude = 20*np.log10(np.abs(magnitude))
magnitude = np.clip(magnitude, -40, 200)
print("pcolormesh call")
#print(type(magnitude))
#print(magnitude.shape)

#maxx = 0.0
#for row in magnitude:
#    for elem in row:
#        if elem> maxx:
#            maxx = elem
#print(maxx)
#sys.exit(-1)
ax.pcolormesh(time, frequencies, magnitude, vmin=0, vmax=50,cmap="bwr",shading="gouraud")
ax.set_ylabel(r"Frequency / Hz",fontsize=20)
ax.set_xlabel(r"Time / s",fontsize=20)
#ax.set_ylim(60,900)#

plt.tight_layout()
plt.savefig("frequencies.png")
