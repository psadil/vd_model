function[p] = listLength_runSim(p)

p.root = [pwd, '\'];



%% Initialise, pretrain, and save weight matrix
[p,weights] = pretrain(p);


%% create stim sequence
[p, stims] = createListLengthStimuli(p);


%% begin sessions

fprintf ('\nThere will be %d sessions in total.\n', p.nSess);

startTime=GetSecs;
for sess = 1:p.nSess,
    
    
    if sess <= p.nSess/2;
        p.layer = 1;
        p.stimCond = sess;
    else
        p.layer = 2;
        p.stimCond = sess - p.nSess/2;
    end
    p.stimSet = 1;
    
    
    %% load session based variables
    p.nRows = p.numRows;
        
    %number of grids in layer
    p.numGrids=p.nGrids(1:p.layer);
    
    % tally of activation by trial and layer
    p.peak_act = zeros(p.nTrials(p.stimCond),2);
    p.totalAct = zeros(p.nTrials(p.stimCond),2);
    
    p.recognition = zeros(p.nTrials(p.stimCond),1);
    p.recogByLayer = zeros(p.nTrials(p.stimCond),p.layer);
    
    %% say what about to happen
    fprintf('\n\nSESSION %d, RAT %d\n', sess, p.ratNum);
    
    %% execute model code
    p = listLengthModel(p,stims,weights);
    
    
    %% save it up
    
    fName = strcat(p.dataDir, '/Session',num2str(sess),'_Rat',num2str(p.ratNum),'.mat');
    save(fName,'p');
    
end

p.runningTime = GetSecs-startTime;

fprintf ('\n\nFinished. \r');
end