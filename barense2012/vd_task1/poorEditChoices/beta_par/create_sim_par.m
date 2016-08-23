% function dPrimePredictions = create_sim_simplex(firstRat, lastRat, etaExp, G_exp)

function [] = create_sim_par(firstRat, lastRat)



% name of folder containing stimuli to use



% data is in DIFFERENCE of dPrime, from first to second half of trials.
% order of [caudal_LA, caudal_HA, PRC_LA, PRC_HA]
% want data: caudal LA => no diff
%            caudal HA => dPrime at 2, then 0 (diff of 2?)
%            PRC LA => no diff
%            PRC HA => no diff
% data = [0, -1.5, 0, 0];

% initialize starting value of parameters
% etaExp
% G_exp
%%% numTrainCycles
%%% numEncodingCycles
% startParms = [ .01 , .8 ]; % .09
% startParms = [ .7 , 20 ];
%_genLog', '_a', num2str(a), ...
% '_b', num2str(b), '_v', num2str(v)];
%% Creates a run_sim file for a single, yoked-pair simulation.

%% Set all variables to default values, which can be overwritten subsequently.


% save RUN_SIM.mat RUN_SIM;

% create stimuli for use
A = .6; % .7
B = .3; % .6
train = 500;
eta = train^-A;
g = .5+10*train^-B;
k = .08;
noise = .5e-6;
leng = 6;
startCrit = 0;


startParms = [ eta , g, k ];

etaExp = startParms(1);
G_exp = startParms(2);
k_expt = startParms(3);


% exptName = 'october17_2015';
nameOfFolder = ['eta', num2str(eta), '_g', num2str(g), ...
    '_K', num2str(k), '_A', num2str(A) ,'_B', num2str(B), '_20enc20_', num2str(train), ...
    'trnNOrand_','5pk_20Fix_NOupOnSwitch_', num2str(noise),'noisThresh_', num2str(startCrit), 'stCrt_',num2str(leng), 'leng_noABS'];

preInitial;

if exist('p', 'var')
    clear p
end

% parpool('AttachedFiles', 'RUN_SIM.mat');

parfor rat = firstRat:lastRat
    
    p = struct();
    
%     A = .6; % .7
%     B = .3; % .6
%     train = 500;
%     eta = train^-A;
%     g = .5+10*train^-B;
%     k = .08;
%     noise = .5e-6;
%     leng = 6;
%     startCrit = 0;
%     
%     
%     startParms = [ eta , g, k ];
%     
%     etaExp = startParms(1);
%     G_exp = startParms(2);
%     k_expt = startParms(3);
    
    % p.exptName = 'october17_2015';
%     p.nameOfFolder = ['eta', num2str(eta), '_g', num2str(g), ...
%         '_K', num2str(k), '_A', num2str(A) ,'_B', num2str(B), '_20enc20_', num2str(train), ...
%         'trnNOrand_','5pk_20Fix_NOupOnSwitch_', num2str(noise),'noisThresh_', num2str(startCrit), 'stCrt_',num2str(leng), 'leng_noABS'];
    
    
    % even at 200 rows (possible 200^2 unique stimuli), that's not enough
    % to contain the 16^4 possible part combinations
    
    p.numRows = 200; %variables with 'num' to denote number are used to define RUN_SIM matrix (and translated to another name before used in simulation)
    p.numLayers = 2;
    
    p.numGrids_Caudal = 4;
    p.numGrids_PRC = 1;
    p.nGrids = [p.numGrids_Caudal, p.numGrids_PRC];
    p.maxNumGrids = max(p.nGrids);
    p.nStimFactors = 4; % number of levels for each dimension
    
    p.components = 8;
    p.numInputDims_Caudal = p.components/p.numGrids_Caudal;   % Hm...assume that all dimensions are equally sampled across caudal grids...
    p.numInputDims_PRC = p.components;
    p.numInputDims = [p.numInputDims_Caudal, p.numInputDims_PRC];   % should be back to [15 3 0] (using only the first 15 components)
    
    p.decision_noise = noise;
    p.maxFixations = [20, 25]; % should it be based on empirical data? total # saccades on match trials = 20
    % first == low ambig, second == high ambig
    % should be [20 25]
    p.k_expt = k_expt;
    p.A = A; % was 0.8 %% Pre-training parameter. The bigger A is, the faster ETA decreases, and the smaller the amount of learning on the weights for all units.
    %     p.A_encoding = .1; % was 2
    p.etaExp = etaExp;
    p.B = B; %was .8 Pre-training parameter. The bigger B is, the faster G decreases, and the smaller the neighbourhood of the winner that gets updated.
    p.G_exp = G_exp;
    %     p.sigma = .1;  % currently, these don't change anything...
    %     p.eta = .01; % currently, these don't change anything...
    p.numTrainCycles = [train, train];
    p.numEncodingCycles = [20, 20]; % now better described as encoding cycles per fixation [LA, HA]
    p.numFeaturesToSample = [p.numGrids_Caudal,p.numGrids_Caudal]; % first == lesion, second == control
    p.fixn_ratio_lowHigh = [.3, .5]; % now describes ratio of within/total
    p.outsideRatio = [.2,.1];
    p.sizeOfPeak = 5;
    p.filtPeak = p.numRows+1;
    p.fives = 0;
    p.variableEncode = 1;
    p.diffEncode = 1;
    p.numThresh = 2;
    p.lengthOfCrit = leng;
    p.famil_diff_thresh_start=[startCrit; startCrit];
    p.setPre = 0;
    p.nameOfFolder = nameOfFolder;
    
    
    p.totalInpDimsConditions = 2; %%Once with small DIMS for caudal, once each with small and large DIMS for intact
    %% Set the numb er of conditions to be 1 by default, for all conditions except INPUT_DIMS.
    p.totalTrainingConditions = 1; %number of different 'pretraining phase' lengths to run
    p.totalStimConditions = 1;
    p.totalStimSets = 1;
    p.totalSimulations = 1;
    p.totalEncodingConditions = 1;
    %     p.expt = 'null'; %% For when running a single session
    p.nMismatch = 36; % for now, define here the number of each trial type for use in actual simulations (i.e., cannot vary in RUN_SIM matrix)
    p.nMatch = 36;
    
    p.nTrials = p.nMismatch+p.nMatch;
    
    %%Make the grid_matrix for calculating city-block distance later.
    [cols, rows] = meshgrid(1:p.numRows);
    p.gridMat = cat(3, rows, cols);
    %% NEED TO SPECIFY WHETHER POSTERIOR OR PRC, so can grab correct element of thresh_start
    
    p.famil_diff_thresh=[repmat(p.famil_diff_thresh_start(1),1,p.lengthOfCrit); repmat(p.famil_diff_thresh_start(2),1,p.lengthOfCrit)];
    % p.famil_diff_thresh2=ones(1,3)*p.famil_diff_thresh_start(2);
    % p.famil_diff_thresh=[p.famil_diff_thresh1, p.famil_diff_thresh2];
    % p.famil_diff_thresh=p.famil_diff_thresh(randperm(length(p.famil_diff_thresh)));
    
    
    p.comparedFeat = zeros(p.nTrials, p.numGrids_Caudal);
    
    fprintf('Creating a new experiment...\n');
    p.expt = 'october17_2015'; %input('\nEnter experiment name: ', 's');
    p.totalStimConditions = 2;%input('\nEnter no. of stimulus conditions: ');
    p.totalStimSets = 1;%input('\nEnter no. of stimulus sets: ');
    
    
    test = load('RUN_SIM.mat');
    
    RUN_SIM = test.RUN_SIM_pre;
    
    run_sim_par(rat,RUN_SIM,p)
    
end

plotFamilDiffs(1, lastRat, nameOfFolder);
