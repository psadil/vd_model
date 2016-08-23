% function dPrimePredictions = create_sim_simplex(firstRat, lastRat, etaExp, G_exp)

function [] = create_sim_simplex(firstRat, lastRat)

global consts


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
A = .6; % .7
B = .3; % .6
train = 500;
eta = train^-A;
g = .5+10*train^-B;
k = .08;
noise = 0;
leng = 6;
startCrit = .5e-5;

% testing out with binomail distribution
prob = .25;

%% parameters for generalized logistic function activation
% 
% a = 0; % lower bound
% k = 1; % upper bound
% c = 1; % typically 1
% q = 10; % related to y(0)
% b = 1; % growth rate
% v = 1; % place of maximal growth

%%


startParms = [ eta , g, k ];

etaExp = startParms(1);
G_exp = startParms(2);
k_expt = startParms(3);

consts.exptName = '16dec2015';
consts.nameOfFolder = ['eta', num2str(eta), '_g', num2str(g), ...
    '_K', num2str(k), '_A', num2str(A) ,'_B', num2str(B), '_20enc20_', num2str(train), ...
    'trn_','5pk_20Fix_', num2str(noise),'noisThresh_', num2str(startCrit), 'stCrt_',num2str(leng), 'leng_', ...
    num2str(prob), 'prob']; %_genLog', '_a', num2str(a), ...
   % '_b', num2str(b), '_v', num2str(v)];
%% Creates a run_sim file for a single, yoked-pair simulation.

%% Set all variables to default values, which can be overwritten subsequently.


for rat = firstRat:lastRat
    
    if exist('p', 'var')
        clear p
    end
    
%     
%     p.a = a; % lower bound
%     p.k = k; % upper bound
%     p.c = c; % typically 1
%     p.q = q; % related to y(0)
%     p.b = b; % growth rate
%     p.v = v; % place of maximal growth
%     
    
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
    p.numEncodingCycles = [80, 80]; % now better described as encoding cycles per fixation [LA, HA]
    p.prob = prob;
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
    p.nameOfFolder = consts.nameOfFolder;
    
       
    p.totalInpDimsConditions = 2; %%Once with small DIMS for caudal, once each with small and large DIMS for intact
    %% Set the numb er of conditions to be 1 by default, for all conditions except INPUT_DIMS.
    p.totalTrainingConditions = 1; %number of different 'pretraining phase' lengths to run
    p.totalStimConditions = 1;
    p.totalStimSets = 1;
    p.totalSimulations = 1;
    p.totalEncodingConditions = 1;
    p.expt = 'null'; %% For when running a single session
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
    p.expt = consts.exptName; %input('\nEnter experiment name: ', 's');
    p.totalStimConditions = 2;%input('\nEnter no. of stimulus conditions: ');
    p.totalStimSets = 1;%input('\nEnter no. of stimulus sets: ');
     
    
    %% Some numbers I'll need in this function
    %PUT IN p IF IT TURNS OUT I NEED THEM IN OTHER FUNCTIONS
    numTrials = p.totalStimConditions*p.totalStimSets; %% i.e. number of trials per "rat"
    totalRepeats = p.totalStimConditions*p.totalTrainingConditions;
    
    %% Save 'p' structure for run_sim.m and write_spreadsheet.m functions
    save SIM_PARAMS.mat p
        
    %%%%%%%%%% Now set up the RUN_SIM matrix %%%%%%%%%%%%
    %% Create EXPT column of RUN_SIM array
    expt(2:numTrials+1,1) = {p.expt};
    expt(1,1) = {'EXPT'};
    
    %% Create STIM_CONDITION column of RUN_SIM array
    for fill_up = 1: p.totalStimConditions,
        stim_condition((fill_up-1)*p.totalStimSets+2:fill_up*p.totalStimSets+1, 1) = {fill_up}; %% plus 2 to leave first element for col heading
    end
    
    stim_condition(1,1) = {'STIM_COND'};
    
    %% Create STIM_SET column of RUN_SIM array
    for fill_up_again = 1: p.totalStimConditions,
        for fill_up = 1: p.totalStimSets,
            stim_set((fill_up_again-1)*p.totalStimSets+fill_up+1,1) = {fill_up};
        end
    end
    
    stim_set(1,1) = {'STIM_SET'};
    
    %% Create next column of RUN_SIM array
    for fill_up = 1: numTrials,
        num_rows(fill_up+1,1) = {p.numRows};
    end
    num_rows(1,1) = {'NUM_ROWS'};
    
    %% Squidge together the columns we have so far...
    interim_matrix_1 = cat(2, expt, stim_condition, stim_set, num_rows);
    
    %% And multiply the matrix by the number of other conditions we have
    [r, c] = size(interim_matrix_1);
    r = r-1;
    
    %% Multiply matrix by totalTrainingConditions
    interim_matrix_2 = repmat(interim_matrix_1(2:end,:), [p.totalTrainingConditions 1]);
    interim_matrix_2 = cat(1,interim_matrix_1(1,:),interim_matrix_2); %put col headers back on
    %% And make a new column
    for trainCond = 1: p.totalTrainingConditions,
        num_training_cycles((trainCond-1)*r+1+1:trainCond*r+1,1) = {p.numTrainCycles(trainCond)};
    end
    interim_matrix_2 = cat(2, interim_matrix_2, num_training_cycles);
    
    [r, c] = size(interim_matrix_2);
    r = r-1;
    
    %% Multiply matrix by totalEncodingConditions
    interim_matrix_3 = repmat(interim_matrix_2(2:end,:), [p.totalEncodingConditions 1]);
    interim_matrix_3 = cat(1,interim_matrix_2(1,:),interim_matrix_3); %put col headers back on
    %% And make a new column
    for encodCond = 1: p.totalEncodingConditions,
        num_encoding_cycles((encodCond-1)*r+1+1:encodCond*r+1,1) = {p.numEncodingCycles(encodCond)};
    end
    interim_matrix_3 = cat(2, interim_matrix_3, num_encoding_cycles);
    
    [r, c] = size(interim_matrix_3);
    r = r-1;  % I don't understand the use of these row/column variables, let alone their reduction...
    
    
    n = length(interim_matrix_1(1,:));
    
    %% Multiply matrix by TOTAL_INPUT_DIMS_CONDITIONS
    interim_matrix_4 = repmat(interim_matrix_3(2:end,:), [p.totalInpDimsConditions 1]);
    interim_matrix_4 = cat(1,interim_matrix_3(1,:),interim_matrix_4); %put col headers back on
    
    for InpDims = 1: p.totalInpDimsConditions,
        %% And make a new column for p.numInputDims
        num_input_dims((InpDims-1)*r+1+1:InpDims*r+1,1) = {p.numInputDims(InpDims)};
        %% And a new column to say which "group and layer" of the whole simulation this session is running
        which_gp_lay((InpDims-1)*r+1+1:InpDims*r+1,1) = {InpDims};
    end
    interim_matrix_4 = cat(2, interim_matrix_4, num_input_dims, which_gp_lay);
    
    %% Now add new headings to columns without them.
    interim_matrix_4(1, 1:n) = interim_matrix_1(1,1:n);
    [r, c] = size(interim_matrix_4);
    r = r-1;
    interim_matrix_4{1,n+1} = 'NUM_TRAIN_CYC';
    interim_matrix_4{1,n+2} = 'NUM_ENCOD_CYC';
    interim_matrix_4{1,n+3} = 'NUM_INP_DIMS';
    interim_matrix_4{1,n+4} = 'WHICH_GP_LAYER';
    RUN_SIM = interim_matrix_4;
    
    
    save RUN_SIM.mat RUN_SIM;
    
    %% create stimuli for use
    
    [p,~] = VDcreateStimuli_forC(p);
    
    %
    run_sim(rat)
    
end

plotFamilDiffs(1, lastRat, consts.nameOfFolder);
