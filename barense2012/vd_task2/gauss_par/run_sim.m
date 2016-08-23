function run_sim(p)

% run_sim - runs the simulation defined in the RUN_SIM.mat file

p.root = [pwd, '\'];

dataDir = strcat(p.root,'rats');
if ~exist(dataDir, 'dir'),
    mkdir(dataDir);
end

dataDir = strcat(p.root,'rats/rat',num2str(p.ratNum));
if ~exist(dataDir, 'dir'),
    mkdir(dataDir);
end

% unique stims
[p, stims] = VDcreateStimuli_simple(p);


%% begin sessions

fprintf ('\nThere will be %d sessions in total.\n', p.nSess);

% preTrain this rat, for use in all conds
[p,weights_pre] = VD_pretrain(p);

startTime=GetSecs;
for sess = 1:p.nSess,
    
    
    %% load session based variables
    p.stimCond = mod(sess,3);
    if p.stimCond == 0
        p.stimCond = 1;
    end
    
    % create stim sequence
    [p] = VD_createStimOrder(p);
    
    %% load session based variables
    
    % block HA and LA conditions, useful for when not resetting weights
    p.sess = sess;
    weights = weights_pre;
    if sess > 3
        if sess == 4
            weights = weights_pre;
%         else
            
        end
        p.which_gp_layer = 2 ;
    else
        if sess == 1
            weights = weights_pre;
        end
        p.which_gp_layer = 1;
    end
    
    
    p.stimSet = 1;
    p.nRows = p.numRows;
    p.nTrainCycles = p.numTrainCycles;
    p.numLayers = p.which_gp_layer; % whether PRC is available
    p.nInpDims = p.numInputDims(p.which_gp_layer);
    
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
    
    % last two are: row/col and prev/new
    p.winning = zeros(p.layer,max(p.numGrids),p.nTrials,2,2);
    
    p.nEncodCycles = p.numEncodingCycles(p.stimCond);
    
    %% say what about to happen
    fprintf('\n\nSESSION %d, RAT %d\n', sess, p.ratNum);
    
    %reset famil_diff_thresh so that it's not carried over from previous session
    p.famil_diff_thresh=[repmat(p.famil_diff_thresh_start(1),1,p.lengthOfCrit); repmat(p.famil_diff_thresh_start(2),1,p.lengthOfCrit)];
    
    %% execute model code
    [p, ~] = visDiscrimModel(p,stims,weights);
    
    
    %% save it up
    outerDir = strcat(pwd, '/graphsAndSession');
    if ~exist(outerDir, 'dir'),
        mkdir(outerDir);
    end
    
    dataDir = strcat(outerDir,'/',p.nameOfFolder);
    if ~exist(dataDir, 'dir'),
        mkdir(dataDir);
    end
    
    
    fName = strcat(dataDir, '/Session',num2str(sess),'_Rat',num2str(p.ratNum),'.mat');
    save(fName,'p');
    
    
end

p.runningTime = GetSecs-startTime;

fprintf ('\n\nFinished. \r');

end