function p = tUniqueModel(p,stims,weights_before)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code for making a Self Organising Feature Map grid (SOFM)
%
% Rosie Cowell (Dec 2011)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load the pretrained weights at beginning of visual discrimination task
weights = weights_before;

%% gather stims together
if p.stimCond == 1
    stimuli = stims.LA_misMatch;
elseif p.stimCond == 2
    stimuli = stims.LA_match;
end

%% begin trial loop

for trial = 1:(p.nTrials/2),
    
    %----------------------------------------------------------------------
    %%% Get the two stimuli for this simultaneous visual discrimination trial
    stimPair = squeeze(stimuli(trial,:,:));

            
    % flag to say: update weights
    forTrain = 1;
    % present initial sample stimulus
    [weights, ~, p] = ...
        present_stimulus(stimPair(:,1), weights, p, trial, forTrain);
    
    % simulate delay for appropriate cycles (always 200)
    [p, weights] = interfere(p, weights);
    
    forTrain = 0; % allow stimuli to accumulate interference
    selec_forComp = zeros(p.numLayers,max(p.nGrids),2);
    % re-present initial sample stimulus
    [~, selec_forComp(:,:,1), p] = ...
        present_stimulus(stimPair(:,1), weights, p, trial, forTrain);
    
    % present initial novel stimulus
    [~, selec_forComp(:,:,2), p] = ...
        present_stimulus(stimPair(:,2), weights, p, trial, forTrain);
    
    % calc recognition score
    [p] = calc_recognition(p, selec_forComp, trial);
    
        
end

end
