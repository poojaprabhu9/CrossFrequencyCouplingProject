 function [mival, MImat] = ozkurt_measure(phase_sig, amp_sig, avg)
 %normalised mean vector length, Ozkurt 2010
% function mival = mi_measure(phase_sig, amp_sig)
%
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
    N=length(amp_sig(:,count));
    %Create composite signal
    z = amp_sig(:,count).*exp(1i*phase_sig(:,count));
    MImat(count,1) = (1./sqrt(N)) * abs(mean(z)) / sqrt(mean(amp_sig(:,count).*amp_sig(:,count))); % Normalise
end

if strcmp(avg, 'y')    
   mival=mean(MImat);
    
else
   mival=MImat;
end
