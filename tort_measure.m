function [tortval, MImat] = tort_measure(ph_sig, amp_sig, avg, n_bins)
% Kullback–Leibler Modulation Index, Tort et al., 2008 and 2010
% function tortval = tort_measure(ph_sig, amp_sig, avg, n_bins)
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
% avg - string, either 'y' or 'n', determines whether values are
% averaged over trials or returned as a vector


% % n_bins = 18; % number of phase bins
position=zeros(1,n_bins); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/n_bins;
for j=1:n_bins 
    position(j) = -pi+(j-1)*winsize; 
end

num_trials = size(ph_sig, 2);

parfor i=1:num_trials
    % now we compute the mean amplitude in each phase:
    MeanAmp=zeros(1,n_bins); 
    for j=1:n_bins   
    I = find(ph_sig(:,i) <  position(j)+winsize & ph_sig(:,i) >=  position(j));
    Amp=amp_sig(:,i);
    MeanAmp(j)=mean(Amp(I)); 
    end
     
    % the center of each bin (for plotting purposes) is position+winsize/2
     
    % quantifying the amount of amp modulation by means of a
    % normalized entropy index (Tort et al PNAS 2008):
    
    MI=(log(n_bins)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(n_bins);
    MImat(i,1)=MI;
end

if strcmp(avg, 'y')
   tortval=mean(MImat);
else
   tortval=MImat;
end
