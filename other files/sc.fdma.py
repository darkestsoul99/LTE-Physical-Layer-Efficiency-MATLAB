hit # noinspection PyUnresolvedReferences
import numpy as np
import matplotlib.pyplot as plt

## Functions
def Mod(Numbersymbols,modulation):
    if modulation=="BPSK":
        n=Numbersymbols
        x_bit=np.random.randint(0, 2, n)
        x=np.divide(x_bit*2-1,np.sqrt(2))
    elif modulation=="QPSK":
        n=Numbersymbols*2
        x_bit=np.random.randint(0, 2, n)
        x_aux = np.divide(x_bit*2-1, np.sqrt(2))
        x=x_aux[0:int(x_aux.size/2)]+1j*x_aux[int(x_aux.size/2):]
    elif modulation=="16-QAM":
        n=Numbersymbols*2
        x_bit=np.random.randint(1, 5, n)
        x_aux = np.divide(x_bit*2-5, np.sqrt(10))
        x = x_aux[0:int(x_aux.size / 2)] + 1j * x_aux[int(x_aux.size / 2):]
    elif modulation=="64-QAM":
        n = Numbersymbols * 2
        x_bit=np.random.randint(1, 9, n)
        x_aux = np.divide(x_bit*2-9, np.sqrt(42))
        x = x_aux[0:int(x_aux.size / 2)] + 1j * x_aux[int(x_aux.size / 2):]
    return x, x_bit

def SubMap(x,submap,IFFTlen,FFTlen):
    Q = IFFTlen // FFTlen
    if submap=="Interleaved":
        y = np.zeros(IFFTlen, dtype=complex)
        y[0::Q] = x
    return y

def Upsam(x,os):
    y = np.zeros(len(x)*os, dtype=complex)
    y[0::os] = x
    return y

def RaisedC(ts,Nos,alpha,Trunc):
    Ts = 1
    Nos = os
    T = 1
    t1 = np.arange(-Trunc * Ts, -Ts / Nos + Ts / Nos, Ts / Nos)
    t2 = np.arange(Ts / Nos, Trunc * Ts + Ts / Nos, Ts / Nos)
    t = np.hstack((t1, 0, t2))
    v = np.empty(len(t))
    for i in range(0, len(t)):
        if np.abs(np.abs(alpha * t[i] / T) - 0.5) > 1e-5:
            v[i] = np.sinc(t[i] / T) * np.cos(np.pi * alpha * t[i] / T) / (1 - np.power((2 * alpha * t[i] / T), 2))
        else:
            v[i] = np.sinc(t[i] / T) * np.pi * np.sin(np.pi * alpha * t[i] / T) / (8 * alpha * t[i] / T)
    return v

def Downsam(x,os,SNRdblen):
    y=np.zeros((SNRdblen,x.shape[1]//os),dtype=complex)
    y=x[:,::os]
    return y

def deSubMap(x,submap,IFFTlen,FFTlen):
    y = np.zeros((x.shape[0], FFTlen))
    slice = np.arange(0, IFFTlen, Q)
    y = x[:,slice]
    return y

def DeMod(x,modulation):
    if modulation=="BPSK":
        y = (np.round(np.real(x) * np.sqrt(2)) + 1) // 2
    elif modulation=="QPSK":
        x_aux1 = (np.round(np.real(x)*np.sqrt(2))+1)//2
        x_aux2 = (np.round(np.imag(x)*np.sqrt(2))+1)//2
        y=np.concatenate((x_aux1,x_aux2),1)
    elif modulation=="16-QAM":
        x_aux1 = np.round((np.real(x)*np.sqrt(10)+5)/2)
        x_aux2 = np.round((np.imag(x)*np.sqrt(10)+5)/2)
        y = np.concatenate((x_aux1, x_aux2), 1)
    elif modulation=="64-QAM":
        x_aux1 = np.round((np.real(x)*np.sqrt(42)+9)/2)
        x_aux2 = np.round((np.imag(x)*np.sqrt(42)+9)/2)
        y = np.concatenate((x_aux1, x_aux2), 1)
    return y

## Parameters
#Modulation
mod="BPSK"
#Length of FFT
FFTlen=512
#Number of bits
nbits=FFTlen
#Subcarrier-mapping
submapC="Interleaved"
#SNRdb
SNRdb=np.arange(0,25,1)
#Length of IFFT
IFFTlen=512
#Q
Q=IFFTlen//FFTlen
#Length of CP
CP=20
#Up-sampling
os=4
#Roll-off factor
alpha=0.22
#Truncation
Trunc=8
#Filter
filter=np.array(RaisedC(1,os,alpha,Trunc))
#Number of simulations
Nsim=10**6
Error=np.zeros(SNRdb.size)

for i in range(0,Nsim,1):
    ### System
    ## Transmitter
    #Modulation
    x,x_bit=Mod(nbits,mod)
    #FFT
    FFT_x=np.fft.fft(x,FFTlen)
    #Subcarrier-Mapping
    FFT_inter=SubMap(FFT_x,submapC,IFFTlen,FFTlen)
    #IFFT
    IFFT_x=np.fft.ifft(FFT_inter,IFFTlen)
    #CPIFFTlen
    IFFT_cp=y = np.zeros(IFFTlen+CP, dtype=complex)
    IFFT_cp[0:CP]=IFFT_x[-CP:]
    IFFT_cp[CP:]=IFFT_x
    #Pulse shaping
    #Up-sampling
    IFFT_up=Upsam(IFFT_cp,os)
    #Conv
    y=np.convolve(IFFT_up,filter,'same')
    '''
    plt.plot(IFFT_up)
    plt.plot(y,'r')
    plt.show()
    '''
    ## Channel
    #Noise
    Complex_noise=np.random.rand(len(y))*2-1+np.random.rand(len(y))*2j-1j
    Power_noise=np.power(10,-SNRdb/10)
    r=y+np.sqrt(Power_noise[:,np.newaxis]/Q)*Complex_noise[:,np.newaxis].T
    ## Receiver
    #Down-Sampling
    r_dw=Downsam(r,os,len(SNRdb))
    #Remove CP
    slice=np.arange(0,CP,1)
    r_rcp= np.delete(r_dw,slice,1)
    #Remove-FFT
    r_dFFT=np.fft.fft(r_rcp,IFFTlen,1)
    #Equalization
    r_eq=r_dFFT
    #De-Subcarrier
    r_dsub=np.zeros((SNRdb.size,FFTlen))
    r_dsub=deSubMap(r_eq,submapC,IFFTlen,FFTlen)
    #Remove-IFFT
    r_dIFFT=np.fft.ifft(r_dsub,FFTlen,1)
    #Compute error
    y=DeMod(r_dIFFT,mod)
    X_bit=np.array([x_bit,]*SNRdb.size)
    error=np.sum(y != X_bit,1)
    Error+=error/nbits
Error=Error/Nsim
plt.semilogy(SNRdb,Error)
plt.show()

