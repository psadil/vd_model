function [weights, stopSampling, p, trial_info] = compare_stimuli(stimPair, weights, p, trial_study)

stim = 1;
first_stim_sampled = 1;
stopSampling = 0;
comparison = 0;

%% preparing input stimulus for processing

% stimPair = stimPair(1:p.components,:);
keepSampling = 0;

%%
while ((p.fixations(trial_study) < p.maxFix) && stopSampling == 0) || (keepSampling && stopSampling == 0)  % stop when enough evidence accumulated or maxFix saccades made
    
    comparison = comparison + 1;
    
    %----------------------------------------------------------------------
    % Expose network to this stimulus, calculate selectivity, and update weights
    %----------------------------------------------------------------------
    
    [weights, selec_forComp(:,:,1), p, ~] = ...
                present_stimulus(stimPair(:,stim), weights, p, trial_study);
    
    
    % after completing within-stimulus fixations, some fixations can land,
    % not on the other stimulus, but outside of both
    nothingRatio = rand;
    if nothingRatio > p.outsideRatio(p.stimCond)
        p.fixations(trial_study) = p.fixations(trial_study) + 1; %% total fixations across both stimuli
    end
    
    if first_stim_sampled == 0 
        
        %------------------------------------------------------------------
        % compile structure for determining mismatch
        %------------------------------------------------------------------
        
        judging.selectivity_new = selectivity;
        judging.selectivity_prev = prevStimSelec;
        judging.featuresSampled_new = features_sampled;
        judging.featuresSampled_prev = features_sampled_prev;
        judging.initial_selec = initial_selec;
        
        % storing whether PRC is used
        % second row is always new stim
        p.usePRC(stim,trial_study) = usePRC;
        
        
        %------------------------------------------------------------------
        % determine which layer is most selective
        %------------------------------------------------------------------
        [stopSampling, p, whichCaudal] = determineMisMatch(judging, p, trial_study);
                
        
        %------------------------------------------------------------------
        % save all relevant info. in a structure for criterion analysis
        %------------------------------------------------------------------
        
        trial_info.selectivity=selectivity;
        trial_info.prevStimSelec = prevStimSelec;
        trial_info.initialSelec_prev = prevInitialSelec;
        trial_info.initialSelec_new = initial_selec;
        
        
        
        %------------------------------------------------------------------
        %stop sampling if evidence of mismatch acquired
        %------------------------------------------------------------------
        if stopSampling==1,
            continue;
        end
        
    end
    
    
    %----------------------------------------------------------------------
    %save weights and activations from this stimulus
    %----------------------------------------------------------------------
    
    features_sampled_prev = features_sampled;
    prevStimSelec = selectivity;
    prevInitialSelec = initial_selec;
    
    prevStimActs = acts;
    
    % first row is always previous stim
    p.usePRC(stim,trial_study) = usePRC;
    
    
    %----------------------------------------------------------------------
    %switch to other stim (last thing you do)
    %----------------------------------------------------------------------
    
    stim = (stim+2)/stim - 1; %1 goes to 2, 2 goes to 1
    if first_stim_sampled == 1 && p.fixations(trial_study) >= p.maxFix
        keepSampling = 1;
    elseif first_stim_sampled == 0
        keepSampling = 0;
    end
    first_stim_sampled = 0; %to indicate that we are no longer on the first stim being sampled (so a comparison should be made from now on)
    
end