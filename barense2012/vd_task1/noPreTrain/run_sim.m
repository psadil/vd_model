function run_sim(rat)

% run_sim - runs the simulation defined in the RUN_SIM.mat file
global ROOT;

% ROOT = 'C:/Users/rcowell/Documents/Research/RH_VisDiscrim/';
ROOT = [pwd '\'];

dataDir = strcat(ROOT,'rats');
if ~exist(dataDir, 'dir'),
    mkdir(dataDir);
end

dataDir = strcat(ROOT,'rats/rat',num2str(rat));
if ~exist(dataDir, 'dir'),
    mkdir(dataDir);
end

if nargin == 2
    fprintf('\nScanning file for previous seed.');
    seed_fid = fopen('seed.txt', 'r');
    seed = fscanf(seed_fid, '%f');
else
    fprintf('\nSetting new seed by clock.');
    seed = fix(1e6*sum(clock));
    seed_fid = fopen('seed.txt', 'w');
    fprintf(seed_fid, '%.0f', seed);
    fclose(seed_fid);
end
rand('state', seed);

%create_sim;
%% Don't call the function create_sim.m because have already made RUN_SIM array on local machine
load RUN_SIM.mat RUN_SIM;
load SIM_PARAMS.mat p;



%% Initialise, pretrain, and save weight matrix
[p] = VD_pretrain(p, rat);


%% begin sessions

p.runningTime = zeros(1,rat);
fprintf ('\nThere will be %d sessions in total.\n', size(RUN_SIM,1)-1);
p.ratNum = rat;

startTime=GetSecs;
for sess = 2:size(RUN_SIM, 1),
    
    
    % create stim sequence
    [p] = VD_createStimOrder(p);
    
    %% load session based variables
    
    % block HA and LA conditions, useful for when not resetting weights
    if mod(p.ratNum,2)
        p.stimCond = RUN_SIM{sess,2};
        p.sess = p.stimCond;
    else 
        if any([2,4] == sess)
           p.stimCond = 2;
           p.sess = 1;
        else
            p.stimCond = 1;
            p.sess = 2;
        end
    end
    

    p.stimSet = RUN_SIM{sess,3};
    p.nRows = RUN_SIM{sess,4};
    p.nTrainCycles = RUN_SIM{sess,5};
%     p.nEncodCycles = RUN_SIM{sess,6};
    p.nInpDims = RUN_SIM{sess,7};
    p.which_gp_layer = RUN_SIM{sess,8};
    p.numLayers = p.which_gp_layer; % whether PRC is available
    
    % choose fixation ratio based on stim condition
    p.fixn_ratio = p.fixn_ratio_lowHigh(p.stimCond);

    % simply nicer to refer to things as layers, sometimes
    p.layer = p.which_gp_layer;
    
    % maximum number of fixations allowed (fixations before judgement of
    % 'match' will be made)
    p.maxFix = p.maxFixations(p.stimCond);
    
    % number of features picked up per fixation
    p.nFeaturesToSample = p.numFeaturesToSample(p.layer);
    
    %number of grids in layer
    p.numGrids=p.nGrids(1:p.layer);
    
    % number of input dimensions to those grids
    p.nInpDims=p.numInputDims(1:p.layer);
    
    % variable to decide whether to every use the PRC layer, determined
    % within VD_present_stimulus.m
    p.usePRC = zeros(2,p.nTrials);

    p.nEncodCycles = p.numEncodingCycles(p.stimCond);
    
    %% say what about to happen
    fprintf('\n\nSESSION %d\n', sess-1);
    if RUN_SIM{sess,8} ~= RUN_SIM{sess-1,8} %%If a new network layer
        if p.which_gp_layer == 1 && (p.totalInpDimsConditions == 1 || p.totalInpDimsConditions == 3)
            fprintf('\nGroup "Lesion", Posterior Layer');
        elseif p.which_gp_layer == 1 && p.totalInpDimsConditions == 2
            fprintf('\nGroup "Control", Posterior Layer');
        elseif p.which_gp_layer == 2 && p.totalInpDimsConditions == 3
            fprintf('\nGroup "Control", Posterior Layer');
        elseif p.which_gp_layer == 2 && p.totalInpDimsConditions == 2
            fprintf('\nGroup "Control", Perirhinal Layer');
        elseif p.which_gp_layer == 3
            fprintf('\nGroup "Control", Perirhinal Layer');
        end
        fprintf('\nExpt: %s; StimCondition: %d; StimSet %d; NumInpDims %d; NumRows %d; rat: %d\n', RUN_SIM{sess,1}, RUN_SIM{sess,2}, RUN_SIM{sess,3}, RUN_SIM{sess,7}, RUN_SIM{sess,4}, rat);
    end
    
    %% execute model code
    p = visDiscrimModel(p);
    
    
    %% save it up
    outerDir = strcat(pwd, '/graphsAndSession');
    if ~exist(outerDir, 'dir'),
        mkdir(outerDir);
    end
     
    dataDir = strcat(outerDir,'/',p.nameOfFolder);
    if ~exist(dataDir, 'dir'),
        mkdir(dataDir);
    end
    
    fName = strcat(dataDir, '/Session',num2str(sess-1),'_Rat',num2str(rat),'.mat');
    save(fName,'p');
    
    %reset famil_diff_thresh so that it's not carried over from previous session
    p.famil_diff_thresh=[repmat(p.famil_diff_thresh_start(1),1,p.lengthOfCrit); repmat(p.famil_diff_thresh_start(2),1,p.lengthOfCrit)];
    
end

p.runningTime = GetSecs-startTime;

fprintf ('\n\nFinished. \r');
