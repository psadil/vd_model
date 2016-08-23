function p = listLengthModel(p,stims,weights_before)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code for making a Self Organising Feature Map grid (SOFM)
%
% Rosie Cowell (Dec 2011)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% uncomment when randomizing initial weights
% weights_before = rand(size(weights_before));

%% begin trial loop
for trial_study = 1:p.nTrials(p.stimCond)
    
    stimPair = squeeze(stims(trial_study,:,:,p.stimCond));
    
    % want fresh weights for only the first stim, but updated weights for
    % stims after that
    if trial_study == 1
        % present initial sample stimulus
        [weights, ~, p] = ...
            present_stimulus(stimPair(:,1), weights_before, p, trial_study);
    else
        [weights, ~, p] = ...
            present_stimulus(stimPair(:,1), weights, p, trial_study);
    end
    
end


for trial_test = 1:p.nTrials(p.stimCond)
    
    stimPair = squeeze(stims(trial_test,:,:,p.stimCond));
    
    selec_forComp = zeros(p.numLayers,max(p.nGrids),2);
    % re-present initial sample stimulus
    [~, selec_forComp(:,:,1), p,] = ...
        present_stimulus(stimPair(:,1), weights, p, trial_test);
    
    % present initial novel stimulus
    [~, selec_forComp(:,:,2), p] = ...
        present_stimulus(stimPair(:,2), weights, p, trial_test);
    
    % calc recognition score
    [p] = calc_recognition(p, selec_forComp, trial_test);
    
end

end