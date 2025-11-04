function feature_eeg = EEG_extract_feature(tmp_EEG)
% EEG_extract_feature  Compute basic spectral and entropy-based EEG features.
%   feature_eeg = EEG_extract_feature(tmp_EEG)
%   INPUT:
%       tmp_EEG : EEGLAB-style EEG structure containing:
%           - tmp_EEG.data   : EEG data matrix [channels x time points]
%           - tmp_EEG.srate  : Sampling frequency (Hz)
%           - tmp_EEG.pnts   : Number of time points
%           - tmp_EEG.setname: Dataset name
%   OUTPUT:
%       feature_eeg : Struct containing extracted EEG features
%                     (bandpower, relative power, entropy, kurtosis, etc.)

%   The script divides the EEG signal into 10-second windows (non-overlapping)
%   and extracts spectral power (absolute and relative), Shannon wavelet
%   entropy, kurtosis, and skewness per window.
%
%   Author: Sneha Ray, UCSF (sneharay@ucsf.edu)
%% Initialize feature structure
Fs = tmp_EEG.srate;                     % Sampling frequency
eeg_features = struct('org_set', tmp_EEG.setname);
% Analysis window (in seconds)
time_window = 10;                       
fprintf('Running EEG feature extraction...\n')
%% Sliding window setup
n = 1;
start_idx = 1;
end_idx = start_idx + time_window * Fs;
%% Main feature extraction loop
while true
    % Break if we reach the end of data
    if end_idx > tmp_EEG.pnts
        break
    end
    % Extract segment [channels x samples]
    eeg = tmp_EEG.data(:, start_idx:end_idx);
    % Skip windows with NaNs
    if all(~isnan(eeg(:)))
       %% ---------------- POWER SPECTRAL MEASURES ----------------
        % Compute absolute bandpower for standard EEG bands
        deltapow = bandpower(eeg', Fs, [1 4]);
        thetapow = bandpower(eeg', Fs, [4 8]);
        alphapow = bandpower(eeg', Fs, [8 12]);
        beta1pow = bandpower(eeg', Fs, [12 20]);
        beta2pow = bandpower(eeg', Fs, [20 30]);

        % Total power for normalization (1–30 Hz)
        totalpow = bandpower(eeg', Fs, [1 30]);

        % Relative power (ratio to total power)
        eeg_features.delta_relative(n,:) = deltapow ./ totalpow;
        eeg_features.theta_relative(n,:) = thetapow ./ totalpow;
        eeg_features.alpha_relative(n,:) = alphapow ./ totalpow;
        eeg_features.beta1_relative(n,:) = beta1pow ./ totalpow;
        eeg_features.beta2_relative(n,:) = beta2pow ./ totalpow;

        % Absolute powers
        eeg_features.deltapow(n,:) = deltapow;
        eeg_features.thetapow(n,:) = thetapow;
        eeg_features.alphapow(n,:) = alphapow;
        eeg_features.beta1pow(n,:)  = beta1pow;
        eeg_features.beta2pow(n,:)  = beta2pow;

        %% ---------------- ENTROPY MEASURES ----------------
        % Mean Shannon wavelet entropy across channels
        eeg_features.shannon_entropy(n) = ...
            mean(abs(wentropy(eeg', 'shannon')));
        % Log energy entropy for each channel
        log_entropy = zeros(size(eeg,1), 1);
        for si = 1:size(eeg,1)
            log_entropy(si) = abs(wentropy(eeg(si,:)', 'log energy'));
        end
        eeg_features.log_shannon_entropy(n,:) = log_entropy;
        % Kurtosis per channel
        eeg_features.kurtosis(n,:) = kurtosis(eeg');     
        % Overall skewness of absolute amplitude
        eeg_features.overall_skew_amp(n) = ...
            skewness(reshape(abs(eeg'), 1, []), 0);
    end
    %% Update window indices
    n = n + 1;
    start_idx = end_idx;
    end_idx = start_idx + time_window * Fs;
end
%% Output
feature_eeg = eeg_features;
fprintf('EEG feature extraction completed.\n')
end
