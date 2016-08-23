function [dPrimePredictions] = create_sim_simplex(firstRat, lastRat, parms, nOfFolder)

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
    
    A = .3;
    B = .2;
    train = 0;
    eta = p1;
    %     eta = train^-A;
    sigma2 = 1;
    g = sigma2;
    %     g = .5+10*train^-B;
    k = .25;
    leng = 4;
    startCrit = eta/100;
    %     startCrit = 1e-6;
    noise = startCrit*p2;
    
    
    startParms = [ eta , g, k ];
    
    etaExp = startParms(1);
    G_exp = startParms(2);
    k_expt = startParms(3);
    
    p.exptName = '20feb2016';
    
    p.nameOfFolder = nOfFolder;
    %     p.nameOfFolder = ['eta', num2str(eta), '_g', num2str(g), ...
    %         '_K', num2str(k), '_A', num2str(A) ,'_B', num2str(B), '_20enc20_', ...
    %         '5pk_20Fix_', num2str(noise),'nois_', num2str(startCrit), 'stCrt_',num2str(leng), '_0reload1',...
    %         'altTanh'];
    
    p.nSess = 4;
    p.sigma2 = sigma2;
    p.numRows = 200; %variables with 'num' to denote number are used to define RUN_SIM matrix (and translated to another name before used in simulation)
    p.numLayers = 2;
    
    p.numGrids_Caudal = 4;
    p.numGrids_PRC = 1;
    p.nGrids = [p.numGrids_Caudal, p.numGrids_PRC];
    p.maxNumGrids = max(p.nGrids);
    p.nStimFactors = 4; % number of levels for each dimension
    
    p.components = 8;
    p.numInputDims_Caudal = p.components/p.numGrids_Caudal;
    p.numInputDims_PRC = p.components;
    p.numInputDims = [p.numInputDims_Caudal,  p.numInputDims_PRC];
    p.nDimReps=1; % number of times to repeat dims (cowell2006 == 3)
    
    p.decision_noise = noise;
    p.maxFixations = [20, 25]; % total # saccades on match trials = 20
    p.k_expt = k_expt; % sigmoidal rate param
    p.A = A; % Pre-training parameter. The bigger A is, the faster ETA decreases, and the smaller the amount of learning on the weights for all units.
%         p.a = 1.7159; % tanh param
        p.b = 2/3; % also tanh
    p.a = 150;
    p.b = atanh(2/3); % with this, max value of act will be (2/3)*a=10
    p.etaExp = etaExp;
    p.B = B; % Pre-training parameter. The bigger B is, the faster G decreases, and the smaller the neighbourhood of the winner that gets updated.
    p.G_exp = G_exp;
    p.numTrainCycles = [train, train];
    p.numEncodingCycles = [5, 5]; % now better described as encoding cycles per fixation [LA, HA]
    p.numFeaturesToSample = [p.numGrids_Caudal,p.numGrids_Caudal]; % first == lesion, second == control
    p.fixn_ratio_lowHigh = [.3, .5]; % now describes ratio of within/total
    p.outsideRatio = [.2,.1]; % chance to fixate outside of either stim
    p.sizeOfPeak = 5;
    p.variableEncode = 1;
    p.diffEncode = 1;
    p.numThresh = 2;
    p.lengthOfCrit = leng;
    p.famil_diff_thresh_start=[startCrit; startCrit];
    p.nameOfFolder = p.nameOfFolder;
    
    
    p.totalInpDimsConditions = 2;
    p.totalTrainingConditions = 1; %number of different 'pretraining phase' lengths to run
    p.totalStimConditions = 1;
    p.totalStimSets = 1;
    p.totalSimulations = 1;
    p.totalEncodingConditions = 1;
    p.expt = 'null'; %% For when running a single session
    
    % number of each trial type for use in actual simulations
    p.nMismatch = 36;
    p.nMatch = 36;
    p.nTrials = p.nMismatch+p.nMatch;
    
    %%Make the grid_matrix for calculating city-block distance later.
    [cols, rows] = meshgrid(1:p.numRows);
    p.gridMat = cat(3, rows, cols);
    
    
    fprintf('Creating a new experiment...\n');
    p.expt = p.exptName;
    p.totalStimConditions = 2;
    p.totalStimSets = 1;
    
    
    %% create stimuli for use
    
    run_sim(p);
    
end

%% looking at network outputs


% uncomment if running simplex
dPrimePredictions = calcDPrime(1,lastRat, nOfFolder);

% uncomment if looking at
% plotFamilDiffs(1, lastRat, nOfFolder,0);

end