# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 16:23:49 2023

@author: Luiz
"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

# Generate unevenly spaced data as an example
# Replace this with your own dataset
t = np.sort(np.random.uniform(0, 10, 100))  # Random time points
data = np.sin(t) + np.random.normal(0, 0.1, t.shape)

# Define the passband period range
min_period = 2.0  # Minimum period
max_period = 10.0  # Maximum period


# Calculate the equivalent cutoff frequencies
lowcut = 1 / max_period
highcut = 1 / min_period

# Sample rate (you can set this based on your data)
fs = 100

# Normalize the cutoff frequencies to Nyquist frequency
nyquist = 0.5 * fs
lowcut_normalized = lowcut / nyquist
highcut_normalized = highcut / nyquist

# Design the filter using scipy.signal
b, a = signal.butter(
    4, [lowcut_normalized, highcut_normalized], 
    btype='band')

# Apply the filter to the data
filtered_data = signal.lfilter(b, a, data)

# Plot the original and filtered data
plt.figure()
plt.plot(t, data, 'b-', label='Original Data')
plt.plot(t, filtered_data, 'g-', linewidth=2, label='Filtered Data')
plt.xlabel('Time')
plt.ylabel('Amplitude')
plt.legend()
plt.grid()
plt.show()