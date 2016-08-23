function[p] = listLength_runSim(p)
% listLength_runSim -- controls flow of sessions. That is, sets necessary
% session specific parameters for whichever rat runSim was called for. As
% such, this function might better be defined as a script.

% called by: createSim
% calls: model

% input:
%   p: experimental structure. defines many things. See 'createSim.m' for
%   full list

% output:
%   NA -- could output p, but not used in the majority of cases



%% set up

% Initialise, pretrain, and save weight matrix
% a rat is born
% this rat will experience all sessions
[p,weights] = pretrain(p);

% generate stimuli to be used for all sessions
stims = createListLengthStimuli(p);


%% begin sessions

fprintf ('\nThere will be %d sessions in total.\n', p.nSess);

startTime=GetSecs;
for sess = 1:p.nSess,
    
    
    % stimCondition defines list length. So, when we're working with rats
    % in the control sessions (when layer == 2), stimCond needs to start
    % over again at 1.
    % p.layer stands in for how many layers are available in this session.
    % So, 1 => lesion, 2 => control.
    if sess <= p.nSess/2;
        p.layer = 1;
        p.stimCond = sess;
    else
        p.layer = 2;
        p.stimCond = sess - p.nSess/2;
    end
    
    
    %% load session based variables
    p.nRows = p.numRows;
    
    %number of grids in layer
    p.numGrids=p.nGrids(1:p.layer);
    
    % tally of activation by trial and layer
    p.peak_act = zeros(p.nTrials(p.stimCond),2);
    p.totalAct = zeros(p.nTrials(p.stimCond),2);
    
    %----------------------------------------------------------------------
    % recognition scores are what we're, ultimately, after.
    %----------------------------------------------------------------------
    p.recognition = zeros(p.nTrials(p.stimCond),p.nStimSets);
    p.recognition_gauss = zeros(size(p.recognition));
    p.recogByLayer = zeros(p.nTrials(p.stimCond),p.layer,p.nStimSets);
    p.recogByLayer_gauss = zeros(size(p.recogByLayer));
    
    fprintf('\n\nSESSION %d, RAT %d\n', sess, p.ratNum);
    
    %% execute model
    p = listLengthModel(p,stims,weights);
    
    
    %% save it up
    p.runningTime = GetSecs-startTime;
    fName = strcat(p.dataDir, '/Session',num2str(sess),'_Rat',num2str(p.ratNum),'.mat');
    save(fName,'p');
    
end

fprintf ('\nFinished. \r');
end