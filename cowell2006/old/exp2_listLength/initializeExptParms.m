function [ p ] = initializeExptParms( expt )
%initializeExptParms -- initializes experimental specific parameters


%% first, initialize parameters shared by every expt
p = struct();
p.expt = expt;

%--------------------------------------------------------------------------
% Archetecture of model
%--------------------------------------------------------------------------

% nodes in a row of grid (total grid is numRows x numRows)
p.numRows = 200;
p.numLayers = 2;

p.numGrids_Caudal = 4;
p.numGrids_PRC = 1;
p.nGrids = [p.numGrids_Caudal, p.numGrids_PRC];

p.components = 8; % num elemental features
p.nStimFactors = 4; % number of levels for each dimension
p.numInputDims_Caudal = p.components/p.numGrids_Caudal;
p.numInputDims_PRC = p.components;
p.numInputDims = [p.numInputDims_Caudal, p.numInputDims_PRC];

%--------------------------------------------------------------------------
% pre-training parms
%--------------------------------------------------------------------------

% cycles to go through pretraining (NOTE: here's a place to fix in
% future projects: don't cycle through a set number of times, cycle
% through until a set error has been reached.
p.numTrainCycles = 500;


% The bigger A is, the faster ETA decreases, which makes for slower
% learning that occurs during pretraining
p.A = .6;

% The bigger B is, the faster G decreases, and the smaller the
% neighborhood of the winner that gets updated.
p.B = .3;


%--------------------------------------------------------------------------
% experimental parms
%--------------------------------------------------------------------------

% width of gaussian neighborhood learning equation
p.G_exp = .5+10*p.numTrainCycles^-p.B;

% rate of sigmoid activation function
p.k_expt = .08;

% learning rate. amount that winning node moves closer to input during
% expt
p.etaExp = p.numTrainCycles^-p.A;

% encoding cycles per presentation of stimulus
p.numEncodingCycles = 20;

% define how many nodes to include in selectivity calculation
p.sizeOfPeak = 5;

% used for later calculation of city-block distance.
[cols, rows] = meshgrid(1:p.numRows);
p.gridMat = cat(3, rows, cols);



%% experiment specific parameters

if expt == 1 % delay
    
    p.eta_int = 0.1;
    
    p.delayCycles = [0,200,400,600,800];
    
    p.nSess = length(p.delayCycles) * 2;
    p.nMismatch = 4;
    p.nMatch = 0;
    p.nStimSets = 1;
    
    
elseif expt == 2 % listLength
    
    % different number of stimuli list lengths
    p.nMismatch = [1,6,12,18];
    p.nTrials = p.nMismatch;
    p.nMatch = 0;
    p.nStimSets = 4;
    p.delayCycles = 0;
    
    p.nSess = length(p.nTrials) * p.numLayers;
    
    
elseif expt == 3 % tUnique
    
    p.nMismatch = 30;
    p.nMatch = 30;
    
    p.delayCycles = 200;
    p.nStimSets = 1;
    
    p.nSess = 4;
    
end


%% General parms again

p.nTrials = p.nMismatch + p.nMatch;


end