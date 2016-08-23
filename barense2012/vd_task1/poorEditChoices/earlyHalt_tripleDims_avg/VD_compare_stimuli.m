function [weights, stopSampling, p] = VD_compare_stimuli(stimPair, weights, p, trial)
%VD_compare_stimuli_controllOfFixation called by visDiscrimModel.m

% want this function to: decide which features of the stimulus to present


% 11/2/14 ps in stick && fixation if statement, added a line to clear
%sample_feat before changing them to 1

% 5/18/15 ps
% commented out all assignments of activations (in an effort to speed up
% processing), since not looking at them directly right now.

% fixations = 0;
stim = 1;
first_stim_sampled = 1;
stopSampling = 0;
comparison = 0;
% p.fixations_total(trial) = 0;

%% preparing input stimulus for processing

stimPair = stimPair(1:p.components,:);
keepSampling = 0;

%%
% keepSampling overrides stopSampling, when stopSampling is true (as
% happens on the super rare occasions when too many fixations were made on
% the first stimulus)
while ((p.fixations(trial) < p.maxFix) && stopSampling == 0) || (keepSampling && stopSampling == 0)  % stop when enough evidence accumulated or maxFix saccades made
    
    comparison = comparison + 1;
    sample_feat = zeros(p.numGrids_Caudal,1);
    avail_feat = randperm(p.numGrids_Caudal);
    sample_feat(avail_feat(1:p.nFeaturesToSample)) = 1; %sample features (at random)
    
    % HUH..,MAYBE THIS NEEDS TO BE COMMENTED OUT: FIXATIONS ARE ONLY
    % ENCODING CYCLES...
    %     p.fixations(trial) = p.fixations(trial) + 1; %% total fixations across both stimuli
    
    if ~first_stim_sampled %force at least one feature in newly sampled stim to
        % be compared with features samped from prev
        % stim.
        while ~any(sample_feat(features_sampled_prev))
            sample_feat = zeros(p.numGrids_Caudal,1);
            avail_feat = randperm(p.numGrids_Caudal);
            sample_feat(avail_feat(1:p.nFeaturesToSample)) = 1; %sample features (at random)
        end
    end
    
%     stick = 1;
    %     while stick == 1,
    %         %------------------------------------------------------------------
    %         % decide whether to stick or switch stimulus (using fixation ratio)
    %         %------------------------------------------------------------------
    %
    %         stick_switch = rand;
    %         %         stickProb = p.fixn_ratio/(p.fixn_ratio+1);
    %         if stick_switch <= p.fixn_ratio
    %             stick = 1;
    %         else %if stick_switch >= stickProb
    %             stick = 0;
    %             continue;
    %         end
    %
    %         %------------------------------------------------------------------
    %         %if stick, sample more features
    %         %------------------------------------------------------------------
    %
    %         if stick %&& (p.fixations(trial) <= p.maxFix),
    %             % sample 3 features on each fixation (note: remove one already-sampled feature, if 1st or 2nd fixation on this stim)
    %             if length(avail_feat) > p.nFeaturesToSample, % make sure that enough features are left to sample from
    %                 avail_feat = avail_feat(2:end);
    %                 %NOTE: all features can only be sampled once (features are being totted up over fixations and sampled at the end)
    %
    %                 % HUH..,MAYBE THIS NEEDS TO BE COMMENTED OUT: FIXATIONS ARE ONLY
    %                 % ENCODING CYCLES...
    %                 %     p.fixations(trial) = p.fixations(trial) + 1; %% total fixations across both stimuli
    %
    %                 avail_feat = avail_feat(randperm(length(avail_feat)));
    %                 sample_feat(avail_feat(1:p.nFeaturesToSample)) = 1; %sample features
    %             else
    %                 stick = 0;
    %             end
    %         end
    %     end
    
    p.sample_feat(trial,:) = sample_feat;
    p.featsSampedByComparison(stim,trial, comparison) = sum(sample_feat);
    
    %     % incriment total fixation counter (used to shrink learning parameter)
    %     p.fixations_total = p.fixations_total+fixations;
    
    % Incorporate all evidence together into stimulus presentation; grab stimulus, sample the specified features, fill remaining inputs with 0.5
    features_sampled = find(sample_feat==1);
    
    %next few lines translate the sampled features into correct input dimensions in the 15 element input vector
    inputs = (features_sampled-1)*p.numInputDims_Caudal+1;
    
    inps_tmp = zeros(p.numInputDims_Caudal-1,length(inputs));
    for inp = 1:(p.numInputDims_Caudal-1)
        inps_tmp(inp,:) = inputs + inp;
    end
    inputs = cat(1,inputs,inps_tmp(:));
    nonInputs = ones(p.components,1);
    nonInputs(inputs)=0;
    stimulus = stimPair(:,stim);
    stimulus(nonInputs==1)=0.5;
    
    %----------------------------------------------------------------------
    % Expose network to this stimulus, calculate selectivity, and update weights
    %----------------------------------------------------------------------
    
    pktot.fin_act_peak = zeros(p.numLayers,max(p.numGrids));
    pktot.fin_act_total = zeros(p.numLayers,max(p.numGrids));
    pktot.init_act_peak = zeros(p.numLayers,max(p.numGrids));
    pktot.init_act_total = zeros(p.numLayers,max(p.numGrids));
    
    [weights, selectivity, initial_selec, p, pktot, initial_weights] = ...
        VD_present_stimulus(stimulus, weights, p, trial, pktot);
    
    
    % after completing within-stimulus fixations, some fixations can land,
    % not on the other stimulus, but outside of both
    nothingRatio = rand;
    if nothingRatio > p.outsideRatio(p.stimCond)
        p.fixations(trial) = p.fixations(trial) + 1; %% total fixations across both stimuli
    end
    
    if first_stim_sampled == 0 %|| (first_stim_sampled == 1 && p.fixations(trial) >= p.maxFix)% to indicate that we have another stim's worth of activations
        %         p.fixations(trial) = fixations;
        
        
        
        %------------------------------------------------------------------
        % compile structure for determining mismatch
        %------------------------------------------------------------------
        
        judging.selectivity_new = selectivity;
        judging.selectivity_prev = prevStimSelec;
        judging.featuresSampled_new = features_sampled;
        judging.featuresSampled_prev = features_sampled_prev;
        %         judging.activations_new = activations;
        %         judging.activations_prev = prevStimActs;
        %         judging.initial_acts = initial_acts;
        judging.initial_selec = initial_selec;
        
        % storing whether PRC is used
        % second row is always new stim
        %         p.usePRC(stim,trial) = usePRC;
        
        
        %------------------------------------------------------------------
        % determine which layer is most selective
        %------------------------------------------------------------------
        [stopSampling, p, whichCaudal] = determineMisMatch(judging, p, trial);
        
        
        %------------------------------------------------------------------
        % look at grids
        %------------------------------------------------------------------
        %         if p.layer == 2 && ~mod(trial,10)
        %             plotTempActs(prevStimActs, initial_acts, p.layer,trial, whichCaudal, p.tType(trial))
        %         end
        
        
        
        
        %------------------------------------------------------------------
        % save all relevant info. in a structure for criterion analysis
        %------------------------------------------------------------------
        
        %         trial_info.famil_difference=famil_difference;
        %         trial_info.mean_act = mean_act;
        %         trial_info.activations = activations;
%         trial_info.selectivity=selectivity;
%         trial_info.prevStimSelec = prevStimSelec;
        %         trial_info.prevStimActs = prevStimActs;
        %         trial_info.initialActs_prev = prevInitialActs;
%         trial_info.initialSelec_prev = prevInitialSelec;
        %         trial_info.initialActs_new = initial_acts;
%         trial_info.initialSelec_new = initial_selec;
        %         trial_info.correlation = correlation;
        %         trial_info.prev_mean_act=prev_mean_act;
        
        % save peak and total activation for later plotting
        [p] = save_peakTotActs(pktot, p, trial, whichCaudal);
        
        
        %------------------------------------------------------------------
        %stop sampling if evidence of mismatch acquired
        %------------------------------------------------------------------
        if stopSampling==1,
            weights = initial_weights;
            continue;
        end
        
    end
    
    
    %----------------------------------------------------------------------
    %save weights and activations from this stimulus
    %----------------------------------------------------------------------
    
    features_sampled_prev = features_sampled;
    prevStimSelec = selectivity;
%     prevInitialSelec = initial_selec;
    %     prevInitialActs = initial_acts;
    %     prev_mean_act = mean_act;
    pktot.prevStimFin_act_peak = pktot.fin_act_peak;
    pktot.prevStimFin_act_total = pktot.fin_act_total;
    pktot.prevStimInit_act_peak = pktot.init_act_peak;
    pktot.prevStimInit_act_total = pktot.init_act_total;
    
%     prevStimActs = acts;
    
    % first row is always previous stim
%     p.usePRC(:,trial) = usePRC;
    
    
    %----------------------------------------------------------------------
    %switch to other stim (last thing you do)
    %----------------------------------------------------------------------
    
    stim = (stim+2)/stim - 1; %1 goes to 2, 2 goes to 1
    if first_stim_sampled == 1 && (p.fixations(trial) >= p.maxFix)
        keepSampling = 1;
    elseif first_stim_sampled == 0
        keepSampling = 0;
    end
    first_stim_sampled = 0; %to indicate that we are no longer on the first stim being sampled (so a comparison should be made from now on)
    
end