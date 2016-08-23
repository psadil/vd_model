function p = delayModel(p,stims,weights_before)

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


% Load the pretrained weights at beginning of visual discrimination task


% preTrainedWeights=rand(p.layer,p.nRows,p.nRows,p.numInputDims(p.numLayers),p.numGrids(1));
preTrainedWeights = weights_before;


%% initialize storage variables

% tally of activation by trial and layer
p.peak_act = zeros(p.nTrials,2);
p.totalAct = zeros(p.nTrials,2);

p.recognition = zeros(p.nTrials,1);
p.recogByLayer = zeros(p.nTrials,p.layer);

%% begin trial loop

% for picking out the stim in the stimulus pair...

for trial = 1:p.nTrials,
    
    %----------------------------------------------------------------------
    %%% Get the two stimuli for this simultaneous visual discrimination trial
    stimPair = squeeze(stims(trial,:,:));

        
    pktot.fin_act_peak = zeros(p.numLayers,max(p.numGrids));
    pktot.fin_act_total = zeros(p.numLayers,max(p.numGrids));
    pktot.init_act_peak = zeros(p.numLayers,max(p.numGrids));
    pktot.init_act_total = zeros(p.numLayers,max(p.numGrids));
    
    % flag to say: update weights
    forTrain = 1;
    % present initial sample stimulus
    [weights, ~, p, pktot] = ...
        present_stimulus(stimPair(:,1), preTrainedWeights, p, trial, pktot, forTrain);
    
    % simulate delay for appropriate cycles
    [p, weights] = interfere(p, weights);
    
    forTrain = 0;
    selec_forComp = zeros(p.numLayers,max(p.nGrids),2);
    % re-present initial sample stimulus
    [~, selec_forComp(:,:,1), p, pktot] = ...
        present_stimulus(stimPair(:,1), weights, p, trial, pktot, forTrain);
    
    % present initial novel stimulus
    [~, selec_forComp(:,:,2), p, pktot] = ...
        present_stimulus(stimPair(:,2), weights, p, trial, pktot, forTrain);
    
    % calc recognition score
    [p] = calc_recognition(p, selec_forComp, trial);
    
        
end

end
