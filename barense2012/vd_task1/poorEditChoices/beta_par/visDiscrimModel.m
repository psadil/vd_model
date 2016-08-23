function p = visDiscrimModel(p,w_full,stims)

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

% global ROOT



% Load the pretrained weights at beginning of visual discrimination task

% weights=zeros(p.layer,p.nRows,p.nRows,p.components,p.numInputDims(p.numLayers),p.numGrids(p.layer));
% weights=zeros(p.layer,p.nRows,p.nRows,p.numInputDims(p.numLayers),p.numGrids(1));

weights = w_full;

% for layer=1:p.layer
%     nInpDims=p.numInputDims(layer);
%     for grid=1:p.nGrids(layer)
%         if ~p.setPre
% %             location = strcat(ROOT,'rats/rat', num2str(p.ratNum), '/pretrainedW__layer', num2str(layer), 'grid', num2str(grid),'.mat');
% %             load(location);
%             weights(layer,:,:,1:nInpDims,grid)=w;
%             fclose('all');
%             
%         else
%             weights = rand(size(weights));
%         end
%     end
% end


% Get stimuli from stimulus files in '/p.expt/conditionXX/' directory.
% location = strcat(ROOT, p.expt, '/condition', num2str(p.stimCond),'/stimuli.mat');
% fid = load(location);
% fclose('all');

fid = struct();
if p.stimCond == 1
    fid.stimuli1 = stims.LA_match;
    fid.stimuli2 = stims.LA_misMatch;
elseif p.stimCond == 2
    fid.stimuli1 = stims.HA_match;
    fid.stimuli2 = stims.HA_misMatch;
end


% p.threshForPlotting_PRC = zeros(1,p.nTrials);
% p.activations_PRC = zeros(p.numRows,p.numRows,p.numGrids_PRC,p.nTrials);
% p.threshForPlotting_Caudal = zeros(1,p.nTrials);
% p.activations_Caudal = zeros(p.numRows,p.numRows,p.numGrids_Caudal,p.nTrials);

% p.activations_new = zeros(p.layer,p.numRows,p.numRows,max(p.numGrids),p.nTrials);
% p.activations_prev = zeros(p.layer,p.numRows,p.numRows,max(p.numGrids),p.nTrials);


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

% p.answer = zeros(1,p.nTrials);
% p.correct = zeros(1,p.nTrials);

%% begin trial loop

% for picking out the stim in the stimulus pair...
tTypeCnt = [0 0];

for trial = 1:p.nTrials,
    %%% Determine trial type and increment count
    tType = p.tType(trial);
    tTypeCnt(tType) = tTypeCnt(tType)+1;
    if tType==1,
        tTypeStr='Mismatch';
    elseif tType==2,
        tTypeStr='Match';
    end
%     fprintf('\n\nTrial %d: Stim Cond %d, %s Trial, Rat %d.\n',trial, p.stimCond, tTypeStr, p.ratNum);
    
    %----------------------------------------------------------------------
    %%% Get the two stimuli for this simultaneous visual discrimination trial
    stim_name = sprintf('stimuli%d', tType); %tType==1 is Mismatch, tType==2 is Match
    stimuli = fid.(stim_name);
    stimPair = squeeze(stimuli(p.stimOrder(tTypeCnt(tType),tType),:,:));
    stimPair = stimPair(:,:,:);
    
    %----------------------------------------------------------------------
    %%% Generate series of saccades and present stimuli, for this trial
    [weights, stop_sampling, p, trial_info] = VD_compare_stimuli(stimPair, weights, p, trial); %weights output on previous trial get input on next trial
    
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
    
    
    [p] = determineCriterion(p ,trial);
    
end % End of loop over trials

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



