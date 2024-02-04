function [pacmat, pacmatFull, freqvec_ph, freqvec_amp] = getPACallMeasures (sig_amp, Fs,...
    measure, sig_ph, ph_freq_vec, amp_freq_vec)
%function [pacmat, freqvec_ph, freqvec_amp] = getPACAllMeasures (sig_pac, Fs,...
%   measure, sig_mod, ph_freq_vec, amp_freq_vec)
%
% This function calculates a matrix of PAC values using either the  MI, PV, GLM, KLMI and normalised MVLMI 
% It assumes that the input is a prefiltered signal since this function 
% does not include any filtering. The output is a matrix of PAC values and
% a plot of this matrix.
% This signal can only take single vector signals as input but the
% provision for multiple trials will be added
%
% REQUIRED INPUTS:
%  sig_pac - signal suspected of containing PAC
%  Fs - sampling frequency
%  measure - measure to be used - it should be: 'esc', 'mi' or 'cfc'
% OPTIONAL INPUTS:
%  sig_mod - signal containing modulating frequency; default = sig_pac
%  ph_freq_vec - range of frequencies for modulating signal; default = 1:5:101
%  amp_freq_vec - range of frequencies for PAC signal; default = 1:5:101

% Set up some parameters for clarity/storage
xbins = ceil((max(ph_freq_vec) - min(ph_freq_vec))/(diff(ph_freq_vec(1:2))));
ybins = ceil((max(amp_freq_vec) - min(amp_freq_vec))/(diff(amp_freq_vec(1:2))));
cent_freq_vec = zeros(xbins,1);
cent_freq_vec2 = zeros(ybins, 1);


if isa(sig_amp, 'cell')
    % Data represents filtered signals
    num_trials = size(sig_amp{1,1},2); 
else
    % Data represents signals 
    num_trials = size(sig_amp,2);
    
end

for i =1:xbins
    upper_bin = min(ph_freq_vec)+i*(diff(ph_freq_vec(1:2)));
    lower_bin = upper_bin-(diff(ph_freq_vec(1:2)));
    cent_freq_vec(i) =  lower_bin + floor((upper_bin- lower_bin)/2);
end


for i =1:ybins
    upper_bin = min(amp_freq_vec)+i*(diff(amp_freq_vec(1:2)));
    lower_bin = upper_bin-(diff(amp_freq_vec(1:2)));
    cent_freq_vec2(i) =  lower_bin + floor((upper_bin- lower_bin)/2);
end


freqvec_amp = cent_freq_vec2;
freqvec_ph = cent_freq_vec;


% Calculate PAC measures 
pacmat = zeros(ybins, xbins);  
pacmatFull = zeros(num_trials, ybins, xbins);
    
    for i = 1:ybins
        parfor j = 1:xbins
            % Calculate matrix of PAC values
            if strcmp(measure, 'mi')
                % Pacmat full of raw mi values
                [pacmat(i,j), pacmatFull(:,i,j)] = mi_measure(sig_ph{1,j}, sig_amp{1,i},'y');
            end

            if strcmp(measure, 'tort')
                n_bins=18;
                [pacmat(i,j), pacmatFull(:,i,j)] = tort_measure(sig_ph{1,j}, ...
                    sig_amp{1,i},'y', n_bins);
            end

            if strcmp(measure, 'plv')
                [pacmat(i,j), pacmatFull(:,i,j)] = plv_measure(sig_ph{1,j}, ...
                    sig_amp{1,i},'y');
            end

            if strcmp(measure, 'ozkurt')
               [pacmat(i,j), pacmatFull(:,i,j)] = ozkurt_measure(sig_ph{1,j}, ...
                    sig_amp{1,i},'y');
            end

            if strcmp(measure, 'glm')
                [pacmat(i,j), pacmatFull(:,i,j)] = glm_measure(sig_ph{1,j}, ...
                    sig_amp{1,i},'y');
            end
            
        end
    end
    

