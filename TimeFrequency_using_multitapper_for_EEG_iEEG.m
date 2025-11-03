%%%%%%%%%%%%%% To analyse and plot time Frequency for a ECoG signal%%%%%%%%%%%%
%%%%%%% Author: Sneha Ray
%%%%%%% Date of Development: 1 July 2025
%%%%%%% Input: EEG/iEEG preprocessed/cleaned data. This pipline compute
%%%%%%%        time-frequency per channel
%%%%%%%%%%%%%%%%%%%


clc
clear
data_path = 'D:\EEG\PreProcessedData';
cd (data_path)
% provide pathe where you want to save the results image
results_save_path = 'D:\EEG\PreProcessedData\multitap_res\ch5';
% channel number to be analyzed
Channel=5;
SUBJlist=dir('patient_*.mat');
%% Specifay multitaper paprameter for time-frequency analysis
      Fs=100; %Sampling Frequency
      frequency_range=[0 35]; %Limit frequencies from 0 to 25 Hz
      taper_params=[3 5];     %Time bandwidth and number of tapers
      window_params=[4 1];    % %Window size is 4s with step size of 1s
      min_nfft=0;             %No minimum nfft
      detrend_opt='constant'  %detrend each window by subtracting the average
      weighting='unity'       %weight each taper at 1
      plot_on=true;           %plot spectrogram
      verbose=true;           %print extra info

      %% Compute the multitaper spectrogram
for ns_i=1:length(SUBJlist)
    clc
    close all
        SUBJname=SUBJlist(ns_i).name;
        prefix_name=SUBJname(9:end-4);
        prefix_title=[SUBJname(9:15) ' Channel - ' num2str(Channel)];
        prefix_save_name = [prefix_name '_timefreq.png'];
        prefix_save_name_res = [prefix_name '_timefreq.mat'];
        %load the data
        load(SUBJname); % load data
        data1=data(Channel,:);
        %%
      figure
      [spect,stimes,sfreqs] = multitaper_spectrogram(data1,Fs,frequency_range, taper_params, window_params, min_nfft, detrend_opt, weighting, plot_on, verbose);
      caxis([-40 40])
      %%%
      ytic_val=[0 4 8 12 20 25 35 40 45];
      ylim([0 45]); 
      yticks(ytic_val)
      yticklabels({'0Hz', '4Hz','8Hz','12Hz','20Hz','25Hz','35Hz','40Hz','45Hz'})
      title(prefix_title)
      xlabel('Time')
      set(gca,'FontSize',14, 'FontName','Calibri', 'FontWeight', 'bold')  
       %%
cd (results_save_path)
saveas(gcf,prefix_save_name);
save(prefix_save_name_res,'sfreqs','spect','stimes');
cd (data_path)
end