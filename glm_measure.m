function [glmval, mival] = glm_measure(phase_sig, amp_sig, avg)
% General linear Model PAC method Penny et al., 2008
num_trials = size(phase_sig, 2);

parfor count = 1:num_trials
     % Building Matrix of regressors. Note : glmfit adds a column of 1s
     X = [cos(phase_sig(:,count)), sin(phase_sig(:,count)) ones(size(phase_sig(:,count)))]; 
    
     [beta,~, stats] = glmfit(X,amp_sig(:,count),'normal','constant','off');      % Fit the GLM

     mival(count,1) = 1- sum(stats.resid.^2)/sum((amp_sig(:,count)-mean(amp_sig(:,count))).^2);  % 1-var(stats.resid)/var(amplitude); % Calculate the explained variance
end

if strcmp(avg, 'y')    
   glmval=mean(mival);
    
else
   glmval=mival;
end
    
