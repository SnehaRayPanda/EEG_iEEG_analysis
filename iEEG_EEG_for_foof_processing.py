# -*- coding: utf-8 -*-
"""
@author: Sneha Ray, UCSF
"""

import os
import numpy as np
import scipy.io
import h5py
from mne.time_frequency import tfr_array_multitaper
from fooof import FOOOF
from scipy.io import savemat

def load_mat_file(filepath):
    try:
        mat = scipy.io.loadmat(filepath)
        data = mat['data']
    except NotImplementedError:
        with h5py.File(filepath, 'r') as f:
            data = np.array(f['data']).T
    return data

def run_multitaper_spectrogram(data, sfreq=100.0, fmin=0.5, fmax=20.0, n_freqs=40, decim=10):
    freqs = np.linspace(fmin, fmax, n_freqs)
    ch_data = data[np.newaxis,:,:]
    power = tfr_array_multitaper(ch_data, sfreq=sfreq, freqs=freqs, n_cycles=2,
                                 output='power', time_bandwidth=4.0, decim=decim, verbose=False)
    return [power[0]], freqs

def run_fooof_on_spectrograms(spectrograms, freqs):
    results = {
        'aperiodic_params': [],
        'r2_scores': [],
        'peak_params': []
    }
    # Define gamma frequency mask to INCLUDE
    freq_low = 2
    freq_high = 30

    # Define a frequency mask to EXCLUDE 60/76/92 Hz (e.g., ±2 Hz window)
    notch_low = 58
    notch_high = 62
    for time_idx, spec in enumerate(spectrograms):
        time_aperiodic, time_r2, time_peaks = [], [], []
        if np.any(np.isnan(spec)) or np.any(np.isinf(spec)):
            print(f"Skipping time {time_idx} due to NaNs/Infs in data")
            time_aperiodic.append([np.nan, np.nan])
            time_r2.append(np.nan)
            time_peaks.append(np.array([]))
        else:
            try:
                # Apply frequency range mask
                freq_mask = (freqs >= freq_low) & (freqs <= freq_high)
                # Mask out 60 Hz ±2 Hz (notch around 60 Hz)
                #notch_mask = (freqs < notch_low) | (freqs > notch_high)
                # Combine both masks
                final_mask = freq_mask #& notch_mask
                freqs_gamma = freqs[final_mask]
                spec_gamma = spec[final_mask]

                #fm = FOOOF(peak_width_limits=[2, 15], max_n_peaks=3, verbose=False)
                fm = FOOOF(peak_width_limits=[2, 15], max_n_peaks=3, min_peak_height=0.1, verbose=False)
                fm.fit(freqs_gamma, spec_gamma)
                time_aperiodic.append(fm.aperiodic_params_)
                time_r2.append(fm.r_squared_)
                time_peaks.append(fm.peak_params_)
            except Exception as e:
                print(f"FOOOF failed on time {time_idx}: {e}")
                time_aperiodic.append([np.nan, np.nan])
                time_r2.append(np.nan)
                time_peaks.append(np.array([]))

        results['aperiodic_params'].append(np.array(time_aperiodic))
        results['r2_scores'].append(np.array(time_r2))
        results['peak_params'].append(time_peaks)

    return results

def save_results_as_mat(results, filename, save_dir):
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, filename)
    savemat(save_path, {
        'aperiodic_params': np.array(results['aperiodic_params'], dtype=object),
        'r2_scores': np.array(results['r2_scores'], dtype=object),
        'peak_params': np.array(results['peak_params'], dtype=object)
    })
    print(f"Saved to {save_path}")

def process_directory(save_dir, specdir, specfilename, windowlength=30):
    from scipy.io import loadmat
    data = loadmat(specdir + specfilename)
    spectrograms = data['spect']
    freqs = data['sfreqs'][0]

    spectavg = np.zeros((spectrograms.shape[0], int(spectrograms.shape[1] / windowlength)))
    for i in range(spectavg.shape[1]):
        window = spectrograms[:, windowlength * i:windowlength * (i + 1)]
        spectavg[:, i] = window.mean(axis=1)

    fooof_result = run_fooof_on_spectrograms(spectavg.T, freqs)
    save_results_as_mat(fooof_result, f"{specfilename}_fooof.mat", save_dir)

# Example usage
save_dir = 'D:\\EEG\\PreProcessedData\\multitap_res\\ch5\\multitap_res\\ch1\\foofres\\'
specdir = 'D:\\EEG\\PreProcessedData\\multitap_res\\ch5\\multitap_res\\ch1\\'

mat_files = [os.path.splitext(f)[0] for f in os.listdir(specdir) if f.endswith('.mat')]

for specfilename in mat_files[0:300]:
    print(f'Processing File {specfilename}')
    try:
        process_directory(save_dir=save_dir, specdir=specdir, specfilename=specfilename, windowlength=60)
    except Exception as e:
        print(f"Failed processing {specfilename}: {e}")