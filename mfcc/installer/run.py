import scipy.io.wavfile
import stft
import matplotlib.pyplot  as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy  as np
import sys,os
import seaborn as sbn
import glob
sbn.set_style("whitegrid")
plt.rcParams["xtick.labelsize"]=20
plt.rcParams["ytick.labelsize"]=20

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
all_songs = glob.glob("input_songs/*")
print(all_songs)

clean_songs = []
for song in all_songs:
    splitter = song.split("/")[1]
    if splitter=="README":
        continue
    else:
        clean_songs.append(song)
for song in clean_songs:
    rate, audData = scipy.io.wavfile.read(song)
    #extract the info
    length = int(len(audData))#/265152#int(len(audData))#audData.shape[0] # this is the number of samples,
    #if you divide length by the rateyou get the length in seconds /rate
    #print(length, length/rate)
    #sys.exit(-1)
    channel1 = audData[:,0][0:int(length/4)]
    print(channel1.shape)
    #scipy.io.wavfile.write("study.wav",rate,channel1)
    print("Total number of elements in the signal %d"%length);
    #convert in double format
    channel1 = np.double(channel1)
    #channel2 = audData[:,1]
    windowSize = 8192 #length of the window to analyse with STFT
    hopSize = 8192  #hopsize between windows

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
    magnitude, frequencies = stft.play(signal, rate , windowSize, hopSize)
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
    #we need basically a dataframe where each column is time
    #then frquencies for each composition
    #and magnitude
    #padd the column to have a full dataframe for each song
    print("plotting frequencies...")
    fig = plt.figure(figsize=[15,10])
    ax = fig.add_subplot(111)
    #ax.imshow(np.abs(magnitude))
    print("Computing magnitude in dB...")
    magnitude = 20*np.log10(np.abs(magnitude))
    magnitude = np.clip(magnitude, -40, 200)
    print("pcolormesh call")
    print(type(magnitude))
    print(magnitude.shape)

    maxx = 0.0
    for row in magnitude:
        for elem in row:
            if elem> maxx:
                maxx = elem
    print(maxx)
    #sys.exit(-1)
    ax.pcolormesh(time, frequencies, magnitude, vmin=0, vmax=50,cmap="bwr",shading="gouraud")
    ax.set_ylabel(r"Frequency / Hz",fontsize=20)
    ax.set_xlabel(r"Time / s",fontsize=20)
    ax.set_ylim(60,900)

    plt.tight_layout()
    filename = song.split("/")[1] +".png"
    plt.savefig(filename)
