function [p, weights] = visDiscrimModel(p,stims,weights)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code for making a Self Organising Feature Map grid (SOFM)
%
% Rosie Cowell (Dec 2011)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Represent the SOFM as a 2-D grid of x,y coordinates
% i.e. 3 dimensions in all: Rows, Cols, Slices
% 3rd-Dimension, slice 1: x-coords
% 3rd-Dimension, slice 2: y-coords





fid1 = struct();
fid2=struct();
if p.stimCond == 1
    fid2.stimuli1 = stims.LA_misMatch;
    fid2.stimuli2 = stims.LA_match;
    
    fid1.stimuli1 = stims.HA_misMatch;
    fid1.stimuli2 = stims.HA_match;
elseif p.stimCond == 2
    fid1.stimuli1 = stims.HA_misMatch;
    fid1.stimuli2 = stims.HA_match;
    
    fid2.stimuli1 = stims.HA_misMatch;
    fid2.stimuli2 = stims.HA_match;
end



%% initialize storage variables

% to be plotted with plotFamilDiffs
p.meanSelectivity_caudal_new = zeros(p.nTrials,min(p.nGrids));
p.meanSelectivity_PRC_new = zeros(p.nTrials,min(p.nGrids));
p.meanSelectivity_caudal_prev = zeros(p.nTrials,min(p.nGrids));
p.meanSelectivity_PRC_prev = zeros(p.nTrials,min(p.nGrids));
p.familDiff_caudal = zeros(p.nTrials,min(p.nGrids));
p.familDiff_PRC = zeros(p.nTrials,min(p.nGrids));


p.prevStimInit_act_peak = zeros(p.layer,p.nTrials);
p.prevStimInit_act_total = zeros(p.layer,p.nTrials);
p.prevStimFin_act_peak = zeros(p.layer,p.nTrials);
p.prevStimFin_act_total = zeros(p.layer,p.nTrials);
p.newStimInit_act_peak = zeros(p.layer,p.nTrials);
p.newStimInit_act_total = zeros(p.layer,p.nTrials);
p.newStimFin_act_peak = zeros(p.layer,p.nTrials);
p.newStimFin_act_total = zeros(p.layer,p.nTrials);
% two layers when PRC is available (control sessions)


% diagnostic for tracking number of features sampled during each
% comparison. The first dim (a 2) refers to stim
p.featsSampedByComparison = zeros(2,p.nTrials,p.maxFixations(p.stimCond));

% which features were sampled during a given trial
p.sample_feat = zeros(p.nTrials,max(p.nGrids));

% incriment total fixation counter (used to shrink learning parameter)
p.fixations = zeros(1,p.nTrials);

% tally of activation by trial and layer
p.peak_act = zeros(p.nTrials,2);
p.totalAct = zeros(p.nTrials,2);


% which features were compared (should always be all 4,in this setup)
p.comparedFeat = zeros(p.nTrials, p.numGrids_Caudal);

%% begin trial loop

% for picking out the stim in the stimulus pair...
tTypeCnt = [0 0];

for trial = 1:p.nTrials,
    %%% Determine trial type and increment count
    tType = p.tType(trial);
    if tType < 3
        tTypeCnt(tType) = tTypeCnt(tType)+1;
    else
        tTypeCnt(tType-2) = tTypeCnt(tType-2)+1;
    end
    
    %----------------------------------------------------------------------
    %%% Get the two stimuli for this simultaneous visual discrimination trial
    if tType == 1 || tType == 2
        stim_name = sprintf('stimuli%d', tType); %tType==1 is Mismatch, tType==2 is Match
        stimuli = fid1.(stim_name);
    else
       stim_name = sprintf('stimuli%d', tType-2);
       stimuli = fid2.(stim_name);
    end
    
    if tType < 3
        stimPair = squeeze(stimuli(p.stimOrder(tTypeCnt(tType),tType),:,:));
    else
        stimPair = squeeze(stimuli(p.stimOrder(tTypeCnt(tType-2),tType-2),:,:));
    end
        
    
    %----------------------------------------------------------------------
    %%% Generate series of saccades and present stimuli, for this trial
    [weights, stop_sampling, p, threshUpdater] = VD_compare_stimuli(stimPair, weights, p, trial); %weights output on previous trial get input on next trial
    
    %     p.activations(:,:,:,trial)=activations;
    
    %----------------------------------------------------------------------
    %%% Record the discrimination choice
    %stop_sampling: 1 = mismatch, 0 = match (note: tType==1 is Mismatch; tType==2 is Match)
    p.answer(trial) = stop_sampling; % answer=1 is mismatch; answer=0 is match
    p.correct(trial) = abs(p.answer(trial)-(p.tType(trial)-1)); % correct=1 is correct; correct=0 is incorrect
    
    
    %----------------------------------------------------------------------
    %%% Analyze the trial information
    
    % need the (1) because double [0 0] occasionally appear.
    %     p.famil_difference(trial)=trial_info.famil_difference(1);
    
    %     p.correlation(trial,:) = trial_info.correlation;
    %     p.selectivity_prev(trial) = trial_info.selectivity_prev;
    %     p.selectivity_new(trial) = trial_info.selectivity_new;
    %     p.activations_new(1:p.layer,:,:,:,trial) = trial_info.activations;
    %     p.activations_prev(1:p.layer,:,:,:,trial) = trial_info.prevStimActs;
    
    
    [p] = determineCriterion(p ,trial, threshUpdater);
    
end % End of loop over trials

% location = strcat(p.root,'rats/rat', num2str(p.ratNum), '/W_stimCond', num2str(p.sess), '_layer', num2str(p.which_gp_layer), '.mat');
% save(location, 'weights');


% will eventually go into funciton that collects performance of all rats
% and plots measures of their performance...
p.Acc_firstHalf = sum(p.correct(1:p.nTrials/2))/(p.nTrials/2);
p.Acc_secondHalf = sum(p.correct(p.nTrials/2+1:end))/(p.nTrials/2);

hitRate_first = sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==1))/sum(p.tType(1:p.nTrials/2)==1);  % a 'yes' (mismatch judgement) on trials that were mismatches (ie, p.tType==1)
FARate_first = sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==2))/sum(p.tType(1:p.nTrials/2)==2); % a 'yes' on matching trials

hitRate_second = sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==1))/sum(p.tType(p.nTrials/2+1:end)==1);  % a 'yes' (mismatch judgement) on trials that were mismatches (ie, p.tType==1)
FARate_second = sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==2))/sum(p.tType(p.nTrials/2+1:end)==2); % a 'yes' on matching trials

%--------------------------------------------------------------------------
% d' was measure of decreased performance for patients (should drastically
% drop during second half of trials on rats lacking PRC layer, but not be
% affected by controls
p.dPrime_first = norminv(hitRate_first)-norminv(FARate_first);
p.dPrime_second = norminv(hitRate_second)-norminv(FARate_second);



end