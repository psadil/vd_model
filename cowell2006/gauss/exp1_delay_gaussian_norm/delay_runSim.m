function[p] = delay_runSim(p)

p.root = [pwd, '\'];


%% Initialise, pretrain, and save weight matrix
[p,weights] = pretrain(p);

% create stim sequence
[p, stims] = createDelayStimuli(p);


%% begin sessions

fprintf ('\nThere will be %d sessions in total.\n', p.nSess);

startTime=GetSecs;
for sess = 1:p.nSess,
    
    
    if sess <= length(p.delayCycles);
        p.nInpDims = p.numGrids_Caudal;
        p.which_gp_layer = 1;
        p.stimCond = sess;
    else
        p.nInpDims = p.numGrids_PRC;
        p.which_gp_layer = 2;
        p.stimCond = sess - length(p.delayCycles);
    end
    p.stimSet = 1;
    
    
    %% load session based variables
    p.nRows = p.numRows;
    p.nTrainCycles = p.numTrainCycles;
    
    p.numLayers = p.which_gp_layer; % whether PRC is available
    
    % simply nicer to refer to things as layers, sometimes
    p.layer = p.which_gp_layer;
    
    % number of features picked up per fixation
    p.nFeaturesToSample = p.numFeaturesToSample(p.layer);
    
    %number of grids in layer
    p.numGrids=p.nGrids(1:p.layer);
    
    % number of input dimensions to those grids
    p.nInpDims=p.numInputDims(1:p.layer);
    
    % variable to decide whether to every use the PRC layer, determined
    % within VD_present_stimulus.m
    p.usePRC = zeros(2,p.nTrials);
    
    
    %% say what about to happen
    fprintf('\n\nSESSION %d, RAT %d\n', sess, p.ratNum);
    
    %% execute model code
    p = delayModel(p,stims,weights);
    
    
    %% save it up
    p.weights = weights;
    fName = strcat(p.dataDir, '/Session',num2str(sess),'_Rat',num2str(p.ratNum),'.mat');
    save(fName,'p');
    p.weights = [];
end

p.runningTime = GetSecs-startTime;

fprintf ('\n\nFinished. \r');
end