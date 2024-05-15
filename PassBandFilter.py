import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import freqz
from scipy.signal import butter, lfilter

def butter_bandpass(lowcut, highcut, fs, order=5):
    return butter(order, [lowcut, highcut], fs=fs, btype='band')

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y





'''
# Plot the frequency response for a few different orders.
plt.figure(1)
plt.clf()
for order in [3, 6, 9]:
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    w, h = freqz(b, a, fs=fs, worN=2000)
    plt.plot(w, abs(h), label="order = %d" % order)

plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],
         '--', label='sqrt(0.5)')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Gain')
plt.grid(True)
plt.legend(loc='best')

'''

'''
A frequência de amostragem , ou taxa de amostragem, é o número de amostras com 
espaçamento igual por unidade de tempo. Por exemplo, se você tem 96 observações 
igualmente espaçadas por dia, então sua taxa de amostragem é 96/dia, ou 
96/24/3600=0,0011 Hz. Hz, que significa por segundo, é amplamente utilizado 
para taxa de amostragem.

A frequência dos seus dados é diferente. Se o sinal for periódico ou 
aproximadamente periódico (digamos, os padrões de consumo de energia 
semelhantes se repetem todos os dias), a frequência é o número de ocorrência 
desse padrão de repetição por unidade de tempo. 
 
De maneira mais geral, se seus dados não parecerem periódicos, seria útil realizar a Análise de Fourier para entender os componentes de frequência em seus dados. 

A Análise de Fourier (Fourier Transform) é a operação para transformar dados no domínio do tempo no domínio da frequência.


'''

# Sample rate and desired cutoff frequencies (in Hz).
fs = 5000.0
lowcut = 500.0
highcut = 1250.0


# Filter a noisy signal.
T = 0.05
nsamples = T * fs


t = np.arange(0, nsamples) / fs
a = 0.02
f0 = 600.0
x = 0.1 * np.sin(2 * np.pi * 1.2 * np.sqrt(t))
x += 0.01 * np.cos(2 * np.pi * 312 * t + 0.1)
x += a * np.cos(2 * np.pi * f0 * t + .11)
x += 0.03 * np.cos(2 * np.pi * 2000 * t)


plt.figure(2)
plt.clf()
plt.plot(t, x, label='Noisy signal')

y = butter_bandpass_filter(x, lowcut, highcut, fs, order=6)


plt.plot(t, y, label='Filtered signal (%g Hz)' % f0)
plt.xlabel('time (seconds)')
plt.hlines([-a, a], 0, T, linestyles='--')
plt.grid(True)
plt.axis('tight')
plt.legend(loc='upper left')

plt.show()

