function [pacmat, pval, freqvec_ph, freqvec_amp] = getPAC (sig_amp, Fs, filterName, measure, ...
    sig_ph, ph_freq_vec, amp_freq_vec, width, n_surr, epochframes, phFiltorder,ampFiltorder)
% This function calculates a matrix of PAC values using either the  MVLMI, PLV, GLM, KLMI and normalised MVLMI 
 
% It uses shuffled datasets to conduct a significance analysis of the PAC
% values found
%
% Basic function call:
% pacmat = getPACallMeasures(sig_pac, Fs, measure)
%
% REQUIRED INPUTS:
%  sig_pac - signal suspected of containing PAC
%  Fs - sampling frequency
%  measure - measure to be used - it should be:  'mi', 'plv', 'glm', 'tort' or 'ozkurt'

%
% OPTIONAL INPUTS:
%  sig_mod - signal containing modulating frequency; default = sig_pac
%  ph_freq_vec - range of frequencies for modulating signal
%  amp_freq_vec - range of frequencies for PAC signal
%  width - width of the wavelet filter; default = 7
%  nfft - the number of points in fft; default = 200
%  num_surr - the number of surrogate data sets to use during significance



if (size(sig_amp,2) ~= size(sig_ph,2))
    sprintf('Error - Signals must have the same number of trials')
    return
end



% Get number of bins as per phase frequencies and amplitude frequencies
xbins = ceil((max(ph_freq_vec) - min(ph_freq_vec))/(diff(ph_freq_vec(1:2))));
ybins = ceil((max(amp_freq_vec) - min(amp_freq_vec))/(diff(amp_freq_vec(1:2))));


% Get the Phase and amplitude time series based on Filtering methods
if (strcmp(filterName, 'morlet'))
            [filt_sig_amp, filt_sig_ph] = filt_signalsWAV(sig_amp, sig_ph, Fs, ...
                                                ph_freq_vec, amp_freq_vec, measure, width);
elseif (strcmp(filterName, 'FIR')) 
            [filt_sig_amp, filt_sig_ph] = filt_signalsFIR(sig_amp, sig_ph, Fs, ...
                                                ph_freq_vec, amp_freq_vec, epochframes,phFiltorder,ampFiltorder);

elseif (strcmp(filterName, 'MP'))
            [filt_sig_amp, filt_sig_ph] = filt_signalsFIR(sig_amp, sig_ph, Fs, ph_freq_vec, amp_freq_vec, epochframes,phFiltorder,ampFiltorder);

end


% Generate surrogates by shuffling the trial numbers
if n_surr ~= 0

[surr_sig_amp, surr_sig_ph] = generate_surrogates(filt_sig_amp, filt_sig_ph,  n_surr) ;

% PAC measures on surrogate data 

[surr_pacmat, surr_pacmatFull, freqvec_ph, freqvec_amp] = getPACallMeasures(surr_sig_amp, Fs, measure, surr_sig_ph, ph_freq_vec, amp_freq_vec);

end


% Different PAC measures on filetered data

[pacmat, pacmatFull, freqvec_ph, freqvec_amp] = getPACallMeasures(filt_sig_amp, Fs, measure, filt_sig_ph, ph_freq_vec, amp_freq_vec);



% permutation test between surrogate and raw  pac values 
nperm=1000; %number of permutation
for i =1:size(pacmatFull,2)
        for j=1:size(pacmatFull,3)
            x=pacmatFull(:,i,j);
            y=surr_pacmatFull(:,i,j);
            [meandiff,Prob,meandiffs] = permTest_mean(x,y,nperm);

            pval(i,j)=Prob;     % Compute p_value  
               
        end
end




