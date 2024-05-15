import numpy as np

def lanczos_bandpass(in_data, in_fc1, in_fc2, in_n, detrend=False, zero_endpts=False):
    # protect input parameters
    data = in_data.copy()
    fc1 = in_fc1[0]
    fc2 = in_fc2[0]
    flt_n = float(in_n[0])
    int_n = int(in_n[0])
    
    NT = len(data)                     # number of elements of data (as integer)
    flt_NT = float(NT)                 # and as floating type
    dtype = data.dtype                 # type code of input data
    
    if data.ndim != 1:
        raise ValueError("Error--data not vector")
    
    if fc2 <= fc1:
        raise ValueError("Error--fc2 must be larger")
    
    if flt_n < (1.3 / (fc2 - fc1)):
        raise ValueError("Error--n too small")
    
    if detrend:
        # remove linear trend from data
        xval = np.arange(NT)
        fit_param = np.polyfit(xval, data, 1)
        slope = fit_param[0]
        yint = fit_param[1]
        trendline1 = slope * xval + yint
        data -= trendline1
    
    if zero_endpts:
        # make the first and last points of the timeseries 0.0
        data_0 = data[0]
        data_ntm1 = data[NT-1]
        data[0] = 0.0
        data[NT-1] = 0.0
    
    k = np.arange((2 * int_n) + 1) - flt_n  # Variable k=[-n,n]
    k[int_n] = np.nan                       # Inititalize k=0
    
    sigma = np.sin(np.pi * k / flt_n) / (np.pi * k / flt_n)
    wbar_k = (
        (np.sin(2. * np.pi * fc2 * k) / (np.pi * k))
              - 
              (np.sin(2. * np.pi * fc1 * k) / (np.pi * k))) * sigma
    wbar_k[int_n] = 2.0 * (fc2 - fc1)
    
    freq = np.full(NT, np.nan)              # initialize freq vector
    NF = int(np.floor(flt_NT / 2.))         # no. of pos. non-0 freq. pts.
    freq[:NF+1] = np.arange(NF+1) / flt_NT  # fill f=0 and f>0 freq.
    if NT % 2 == 0:                        # if NT is even, fill f<0 freq.
        freq[NF+1:] = -freq[1:NF][::-1]
    else:                                  # if NT is odd, fill f<0 freq.
        freq[NF+1:] = -freq[1:NF+1][::-1]
    
    if not np.isfinite(freq).all():        # error check freq fill
        raise ValueError("Error--bad freq fill")
    
    Rbar_n = np.zeros(NT, dtype=dtype)     # init. smoothed freq. response fctn.
    pos_k = k[int_n+1:]                    # Values of k=[1,n]
    pos_wbar_k = wbar_k[int_n+1:]           # Values of wbar_k corresponding to pos_k
