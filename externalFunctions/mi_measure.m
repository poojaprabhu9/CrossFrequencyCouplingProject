 function [mival, val] = mi_measure(phase_sig, amp_sig,avg)
 % Mean-Vector length modulation index, Canolty et al., 2006
% function mival = mi_measure(phase_sig, amp_sig)
%
% Returns a value for the MI measure calculated between two signals.
%
% INPUTS:
%
% phase_sig - the instantaneous phase values for a signal which has been
% filtered for a lower, modulating frequency, passed as a column vector
%
% amp_sig - the amplitude values for a signal which has been filtered for a
% higher, modulated frequency, passed as a column vector 
%
% Author: Angela Onslow, May 2010


num_trials = size(phase_sig, 2);

parfor count = 1:num_trials
    
    %Create composite signal
    z = amp_sig(:,count).*exp(1i*phase_sig(:,count));
    m_raw(count) = mean(z);  %Compute the mean length of composite signal.      

    val(count,1) = abs((m_raw(count)));
end

if strcmp(avg, 'y')    
   mival=mean(val);
    
else
   mival=val;
end

