import numpy as np
import matplotlib.pyplot as plt

# Function to generate artificial signal with noise and periodic components
def generate_signal(
        duration = 10.0, 
        sampling_rate = 100.0 , 
        periods = [1.0, 3.0, 12.0], 
        amplitudes = [1.0, 0.5, 0.2],
        noise_level = 0.3
        ):
    
    '''
    Duration of the signal in seconds
    Sampling rate in Hz
    Periods of the periodic components in s
    Amplitudes of the periodic components
    Noise level
    
    '''
    
    t = np.arange(0, duration, 1/sampling_rate)  # Time vector
    signal = np.zeros_like(t)
    
    for amplitude, period in zip(amplitudes, periods):
        signal += amplitude * np.sin((2 * np.pi * t) / period)
    
    noise = noise_level * np.random.randn(len(t))
    signal += noise
    
    return t, signal




# Generate the signal
t, signal = generate_signal()

# Plot the generated signal
plt.figure(figsize=(12, 6))
plt.plot(t, signal, label='Noisy Signal')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.title('Artificial Signal with Noise and Periodic Components')
plt.legend()
plt.grid(True)
plt.show()
