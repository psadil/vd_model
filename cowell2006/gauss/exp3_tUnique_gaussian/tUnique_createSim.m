function [] = tUnique_createSim(firstRat, lastRat, parms, nOfFolder)
%% Creates a run_sim file for a single, yoked-pair simulation.

if exist('p', 'var')
    clear p
end

outerDir = strcat(pwd, '/graphsAndSession');
if ~exist(outerDir, 'dir'),
    mkdir(outerDir);
end

rng('shuffle');

p1 = parms(1);
p2 = parms(2);
parfor rat = firstRat:lastRat
    
    p = struct();
    p.ratNum = rat;
    
    A = .2;
    B = .2;
    train = 0;
    etaExp = p1;
    sigma2 = .01;
    G_exp = sigma2;
    k_expt = .25;
    p.eta_int = p2;
    
    p.exptName = '12jan2016';
    p.nameOfFolder = nOfFolder;
    
    p.dataDir = strcat(pwd, '/graphsAndSession/', p.nameOfFolder);
    if ~exist(p.dataDir, 'dir'),
        mkdir(p.dataDir);
    end
    
    
    p.delayCycles = 200;
    p.nSess = length(p.delayCycles) * 4;
    
    % even at 200 rows (possible 200^2 unique stimuli), that's not enough
    % to contain the 16^4 possible part combinations
    
    p.numRows = 200; %variables with 'num' to denote number are used to define RUN_SIM matrix (and translated to another name before used in simulation)
    p.numLayers = 2;
    
    p.numGrids_Caudal = 4;
    p.numGrids_PRC = 1;
    p.nGrids = [p.numGrids_Caudal, p.numGrids_PRC];
    p.maxNumGrids = max(p.nGrids);
    p.nStimFactors = 4; % number of levels for each dimension
    
    p.nDimReps = 1;
    p.components = 8;
    p.numInputDims_Caudal = p.components/p.numGrids_Caudal;
    p.numInputDims_PRC = p.components;
    p.numInputDims = [p.numInputDims_Caudal, p.numInputDims_PRC];
    
    p.k_expt = k_expt;
    p.A = A; % was 0.8 %% Pre-training parameter. The bigger A is, the faster ETA decreases, and the smaller the amount of learning on the weights for all units.
    p.etaExp = etaExp;
    p.sigma2 = sigma2;
    %     p.a = 1.7159; % tanh param
    %     p.b = 2/3; % also tanh
    p.a = 150;
    p.b = atanh(2/3);
    p.B = B; %was .8 Pre-training parameter. The bigger B is, the faster G decreases, and the smaller the neighbourhood of the winner that gets updated.
    p.G_exp = G_exp;
    p.numTrainCycles = [train, train];
    p.numEncodingCycles = 5; % now better described as encoding cycles per fixation [LA, HA]
    p.numFeaturesToSample = [p.numGrids_Caudal,p.numGrids_Caudal]; % first == lesion, second == control
    p.sizeOfPeak = 5;
    p.filtPeak = p.numRows+1;
    p.nameOfFolder = p.nameOfFolder;
    
    p.totalInpDimsConditions = 2; %%Once with small DIMS for caudal, once each with small and large DIMS for intact
    
    p.nMismatch = 30;
    p.nMatch = 30;
    
    p.nTrials = p.nMismatch+p.nMatch;
    
    %%Make the grid_matrix for calculating city-block distance later.
    [cols, rows] = meshgrid(1:p.numRows);
    p.gridMat = cat(3, rows, cols);
    fprintf('Creating a new experiment...\n');
    p.expt = p.exptName;
    
    
    %% create stimuli for use
    
    [~] = tUnique_runSim(p)
    
end

plotRecognition(1, lastRat, nOfFolder);
