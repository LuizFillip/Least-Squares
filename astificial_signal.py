# -*- coding: utf-8 -*-
"""
Created on Thu May 16 08:01:11 2024

@author: Luiz
"""

import numpy as np
import matplotlib.pyplot as plt

# Function to generate artificial signal with noise and periodic components
def generate_signal(duration, sampling_rate, frequencies, amplitudes, noise_level):
    t = np.arange(0, duration, 1/sampling_rate)  # Time vector
    signal = np.zeros_like(t)
    
    for amplitude, period in zip(amplitudes, periods):
        signal += amplitude * np.sin((2 * np.pi * t) / period)
    
    noise = noise_level * np.random.randn(len(t))
    signal += noise
    
    return t, signal

# Parameters
duration = 10.0  # Duration of the signal in seconds
sampling_rate = 100.0  # Sampling rate in Hz
periods = [1.0, 3.0, 5.0]  # Frequencies of the periodic components in Hz
amplitudes = [1.0, 0.5, 0.2]  # Amplitudes of the periodic components
noise_level = 0.3  # Noise level

# Generate the signal
t, signal = generate_signal(duration, sampling_rate, frequencies, amplitudes, noise_level)

# Plot the generated signal
plt.figure(figsize=(12, 6))
plt.plot(t, signal, label='Noisy Signal')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.title('Artificial Signal with Noise and Periodic Components')
plt.legend()
plt.grid(True)
plt.show()
