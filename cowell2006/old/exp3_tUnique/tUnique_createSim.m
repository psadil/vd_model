function [] = tUnique_createSim(firstRat, lastRat, nOfFolder)
%%

% inputs:
%   firstRat: int -- First rat to run. Typically 1
%   lastRat: int -- determines how many rats to run. for efficient use of
%      cores available, should be in multiple of 4
%   nOfFolder: char str -- name of folder to store sessions. to be loaded
%      from when plotting recognition scores. saved in dir
%      'graphsAndSession'

% otuput:
%   NA -- however, tUnique_createSim calls 'plotRecognition.m' which plots
%   recognition scores in a few different ways, and saves those figures in
%   nOfFolder.


%% preliminary setup

% clear old variables, generate folder to store sessions
if exist('p', 'var')
    clear p
end

outerDir = strcat(pwd, '/graphsAndSession');
if ~exist(outerDir, 'dir'),
    mkdir(outerDir);
end

rng('shuffle');


%% begin
parfor rat = firstRat:lastRat
    
    p = struct();
    p.ratNum = rat;
    
    A = .6;
    B = .3;
    train = 500;
    etaExp = train^-A;
    G_exp = .5+10*train^-B;
    k_expt = .08;
    p.eta_int = .05;
    
    
    p.nameOfFolder = nOfFolder;
    
    p.dataDir = strcat(pwd, '/graphsAndSession/', p.nameOfFolder);
    if ~exist(p.dataDir, 'dir'),
        mkdir(p.dataDir);
    end
    
    
    p.delayCycles = 200;
    p.nSess = length(p.delayCycles) * 4;
    
    p.numRows = 200;
    p.numLayers = 2;
    
    p.numGrids_Caudal = 4;
    p.numGrids_PRC = 1;
    p.nGrids = [p.numGrids_Caudal, p.numGrids_PRC];
    p.nStimFactors = 4; % number of levels for each dimension
    
    p.components = 8;
    p.numInputDims_Caudal = p.components/p.numGrids_Caudal;
    p.numInputDims_PRC = p.components;
    p.numInputDims = [p.numInputDims_Caudal, p.numInputDims_PRC];
    
    p.k_expt = k_expt;
    p.A = A;
    p.etaExp = etaExp;
    p.B = B;
    p.G_exp = G_exp;
    p.numTrainCycles = train;
    p.numEncodingCycles = 20;
    p.sizeOfPeak = 5;
    
    
    p.nMismatch = 30;
    p.nMatch = 30;
    p.nTrials = p.nMismatch+p.nMatch;
    
    % Make the grid_matrix for calculating city-block distance later.
    [cols, rows] = meshgrid(1:p.numRows);
    p.gridMat = cat(3, rows, cols);
    fprintf('Creating a new experiment...\n');
    
    
    %% create stimuli for use
    
    [~] = tUnique_runSim(p)
    
end

plotRecognition(1, lastRat, nOfFolder);


end