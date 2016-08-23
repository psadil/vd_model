function [] = listLength_createSim(firstRat, lastRat)
% listLength_createSim -- first function to
% NOTE: each iteration of the simulation is refered to as a 'rat.' One rat
% is run for every session, and is run as though it were a lesioned rat
% (which is to say it doesn't have a PRC grid) and a control rat (the rat
% has BOTH a PRC and caudal grid)

% overall, simulation replicates a Kohonen-style, self-organizing map (SOM)


% called by: YOU
% calls: runSim

% input:
%    firstRat: first rat to run (usually 1)
%    lastRat: last rat to run. defines how many rats are being run. Each
%        rat will run for 8 sessions, 4 lesion and 4 control

% output:
%   NA -- however, commented line at bottom (call to 'plotRecognition.m')
%      will generate summary graphs of the simulation


%% prelim

% p is the main storage structure
if exist('p', 'var')
    clear p
end

% for cleanliness, expts get stored in this folder
outerDir = strcat(pwd, '/graphsAndSession');
if ~exist(outerDir, 'dir'),
    mkdir(outerDir);
end


fprintf('\nSetting new seed by clock.\n');
rng('shuffle');


%% begign running rats
for rat = firstRat:lastRat
    
    p = struct();
    p.ratNum = rat;
    
    % nodes in a row of grid (total grid is numRows x numRows)
    p.numRows = 200;
    p.numLayers = 2;
    
    % different number of stimuli list lengths
    p.nMismatch = [1,6,12,18];
    p.nTrials = p.nMismatch;
    p.nStimSets = 4;
    
    p.nSess = length(p.nTrials) * p.numLayers;
    
    p.numGrids_Caudal = 4;
    p.numGrids_PRC = 1;
    p.nGrids = [p.numGrids_Caudal, p.numGrids_PRC];
    p.nStimFactors = 4; % number of levels for each dimension
    
    p.components = 8;
    p.numInputDims_Caudal = p.components/p.numGrids_Caudal;
    p.numInputDims_PRC = p.components;
    p.numInputDims = [p.numInputDims_Caudal, p.numInputDims_PRC];
    
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
    
    % tanh height
    p.a = 1.7159; 
    
    % tanh rate
    p.b = 2/3; 
    
    % width of gaussian as a selectivity filter
    p.sigma2 = .001;
    
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
    
    
    % tweak as necessary to reflect peculiarities in any given simulation
    p.nameOfFolder = ['eta', num2str(p.etaExp), '_g', num2str(p.G_exp), ...
        '_K', num2str(p.k_expt), '_A', num2str(p.A) ,'_B', num2str(p.B),...
        '_20enc20_', num2str(p.numTrainCycles), 'trn_','1stimSets2'];
    
    
    
    p.dataDir = strcat(outerDir, '\' ,p.nameOfFolder);
    if ~exist(p.dataDir, 'dir'),
        mkdir(p.dataDir);
    end
    
    
    fprintf('Creating a new experiment...\n %s \n', p.nameOfFolder);
    
    
    %% enter into runSim, which controls flow of sessions
    
    [~] = listLength_runSim(p);
    
end

% call following to make graphs

% plotFamilDiffs(1, lastRat, p.nameOfFolder);

end