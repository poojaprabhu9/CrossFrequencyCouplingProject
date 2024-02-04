function [plvval, mival] = plv_measure(ph_sig, amp_sig, avg)
% Phase locking value, Cohen 2008

% function plvval = plv_measure(ph_sig, amp_sig, avg)
%
% Returns a value (or vector of values) for the  measure calculated
% between two signals. Signals may contain one of more trials. Multiple
% trials may be averaged so as to return one pac value or a vector of
% pac values calculated for each trial may be returned, depending on the 
% 'avg' argument. Signals should be passed as column vectors, multiple 
% trials stored as multiple columns.
%
% INPUTS:
% ph_sig - signal filtered for a lower, modulating frequency (e.g. theta
% band oscillations)
%
% amp_sig - signal filtered for a higher, modulated frequency (e.g. gamma
% band oscillations)
%
% avg - string, either 'y' or 'n', determines whether PLV values are
% averaged over trials or returned as a vector
%
% Author: Angela Onslow, May 2010


if size(ph_sig, 2) ~= size(amp_sig, 2)
    sprintf('Error - Signals must have the same number of trials')
    return
end
num_trials = size(ph_sig, 2);



parfor count = 1:num_trials
    
        % Apply PLV algorith, from Cohen et al., (2008)
        amp_phase(:,count) = angle(hilbert(detrend(amp_sig(:,count)))); % Phase of amplitude envelope
        mival(count,:) = abs(mean(exp(1i*(ph_sig(:,count)-amp_phase(:,count)))));
end



if strcmp(avg, 'y')    
   plvval=mean(mival);
    
else
   plvval=mival;
end
