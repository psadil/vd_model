function [] = listLength_createSim(firstRat, lastRat)
%% Creates a run_sim file for a single, yoked-pair simulation.

if exist('p', 'var')
    clear p
end


outerDir = strcat(pwd, '/graphsAndSession');
if ~exist(outerDir, 'dir'),
    mkdir(outerDir);
end


fprintf('\nSetting new seed by clock.\n');
rng('shuffle');


parfor rat = firstRat:lastRat
    
    p = struct();
    p.ratNum = rat;
    
    A = .6;
    B = .3;
    train = 500;
    etaExp = train^-A;
    G_exp = .5+10*train^-B;
    k_expt = .08;
    
    p.nameOfFolder = ['eta', num2str(etaExp), '_g', num2str(G_exp), ...
        '_K', num2str(k_expt), '_A', num2str(A) ,'_B', num2str(B), '_20enc20_', num2str(train), ...
        'trn'];
    
    p.dataDir = strcat(pwd, '/graphsAndSession/', p.nameOfFolder);
    if ~exist(p.dataDir, 'dir'),
        mkdir(p.dataDir);
    end
    
    
    p.numRows = 200; %variables with 'num' to denote number are used to define RUN_SIM matrix (and translated to another name before used in simulation)
    p.numLayers = 2;
    
    % different number of stimuli list lengths
    p.nMismatch = [1,6,12,18];
    p.nTrials = p.nMismatch;
    
    p.nSess = length(p.nTrials) * p.numLayers;
    
    p.numGrids_Caudal = 4;
    p.numGrids_PRC = 1;
    p.nGrids = [p.numGrids_Caudal, p.numGrids_PRC];
    p.nStimFactors = 4; % number of levels for each dimension
    
    p.components = 8;
    p.numInputDims_Caudal = p.components/p.numGrids_Caudal;
    p.numInputDims_PRC = p.components;
    p.numInputDims = [p.numInputDims_Caudal, p.numInputDims_PRC];
    
    p.k_expt = k_expt;
    p.A = A; % was 0.8 %% Pre-training parameter. The bigger A is, the faster ETA decreases, and the smaller the amount of learning on the weights for all units.
    p.etaExp = etaExp;
    p.B = B; %was .8 Pre-training parameter. The bigger B is, the faster G decreases, and the smaller the neighbourhood of the winner that gets updated.
    p.a = 1.7159; % tanh param
    p.b = 2/3; % also tanh
    p.sigma2 = .001;
    p.G_exp = G_exp;
    p.numTrainCycles = train;
    p.numEncodingCycles = 20; % now better described as encoding cycles per fixation [LA, HA]
    p.sizeOfPeak = 5;
        
    
    %%Make the grid_matrix for calculating city-block distance later.
    [cols, rows] = meshgrid(1:p.numRows);
    p.gridMat = cat(3, rows, cols);
    
    
    fprintf('Creating a new experiment...%s\n', p.nameOfFolder);
    
    
    %% create stimuli for use
    
    [~] = listLength_runSim(p)
    
end

% plotFamilDiffs(1, lastRat, consts.nameOfFolder);

end