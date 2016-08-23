function[p] = tUnique_runSim(p)

p.root = [pwd, '\'];


%% Initialise, pretrain, and save weight matrix
[p,weights] = pretrain(p);

% create stim sequence
[p, stims] = createTUniqueStimuli(p);


%% begin sessions

fprintf ('\nThere will be %d sessions in total.\n', p.nSess);

startTime=GetSecs;
for sess = 1:p.nSess,
    % session order:
        % 1) lesion, trial-unique
        % 2) lesion, repeated stims
        % 3) control, trial-unique
        % 4) control, repeated stims
    
    
    if sess <= p.nSess/2;
        p.layer = 1;
        p.stimCond = sess;
    else
        p.layer = 2;
        p.stimCond = sess - p.nSess/2;
    end
        
    %% load session based variables
    p.nRows = p.numRows;
        
    % number of grids in layer
    p.numGrids=p.nGrids(1:p.layer);
            
    p.recognition = zeros(p.nTrials/2,1);
    p.recognitionByLayer = zeros(p.nTrials/2,2);
    
    % say what about to happen
    fprintf('\n\nSESSION %d, RAT %d\n', sess, p.ratNum);
    
    %% execute model code
    p = tUniqueModel(p,stims,weights);
    
    
    %% save it up
    
    fName = strcat(p.dataDir, '/Session',num2str(sess),'_Rat',num2str(p.ratNum),'.mat');
    save(fName,'p');
    
end

p.runningTime = GetSecs-startTime;

fprintf ('\n\nFinished. \r');
end