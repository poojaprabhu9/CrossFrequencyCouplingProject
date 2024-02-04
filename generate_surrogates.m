function [surr_amp_sig, surr_ph_sig] = generate_surrogates(sig_amp, sig_ph,  n_surr) 

% Create matrix for amplitude
surr_amp_sig = [];

% Create matrix of zeros for phase
surr_ph_sig = [];

% warning('Use with caution - NEEDS WORK');

% Get random phase and amplitude trials for each surrogate

% Set seed 
rng('default');
rng(1);

%randomise the trial numbers of phase signal
rand_phase_trial = randi([1 size(sig_ph{1,1},2)], 1, n_surr.*1.5); 

% Set seed if specified
rng('default');
rng(2);

%randomise the trial numbers of amplitude signal
rand_amp_trial = randi([1 size(sig_amp{1,1},2)], 1, n_surr.*1.5);

% Remove instances where the same phase and amplitude trial has
% been selected
find_same_trials = find(rand_amp_trial-rand_phase_trial==0);
rand_phase_trial(find_same_trials) = [];
rand_amp_trial(find_same_trials) = [];

% Now select the correct number of surrogates
rand_phase_trial    = rand_phase_trial(1:n_surr);
rand_amp_trial      = rand_amp_trial(1:n_surr);

disp('Computing surrogate data...');
for n_amp=1:size(sig_amp,2)
    for surr = 1:n_surr
        % Split data
        surr_amp_sig{1,n_amp}(:,surr) = sig_amp{1,n_amp}(:,rand_amp_trial(surr));
        
    end
end

for n_ph=1:size(sig_ph,2)
    for surr = 1:n_surr
        % Split data
        surr_ph_sig{1,n_ph}(:,surr) = sig_ph{1,n_ph}(:,rand_phase_trial(surr));
        
    end
end

end