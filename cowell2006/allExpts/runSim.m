function[p] = runSim(p)
% runSim -- controls flow of sessions. That is, sets necessary
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
stims = init_stimuli(p);


%% begin sessions

fprintf ('\nThere will be %d sessions in total.\n', p.nSess);

startTime=GetSecs;
for sess = 1:p.nSess,
    
    
    % stimCond will index:
    %   nDelay cycles in expt 1
    %   listLength in expt 2
    %   whether stimuli are trials unique (stimCond == 1) or not (stimCond == 2)
    
    % layer indicates lesion or control
    %   1 => lesion
    %   2 => control
    if sess <= p.nSess/2;
        p.layer = 1;
        p.stimCond = sess;
    else
        p.layer = 2;
        p.stimCond = sess - p.nSess/2;
    end
    
    
    
    %----------------------------------------------------------------------
    % recognition scores are what we're, ultimately, after.
    %----------------------------------------------------------------------
    p.recognition = zeros(p.nTrials(p.stimCond),p.nStimSets);
    p.recogByLayer = zeros(p.nTrials(p.stimCond),p.layer,p.nStimSets);
    p.corr = zeros(p.nTrials(p.stimCond),p.nStimSets);
    p.corrByLayer = zeros(p.nTrials(p.stimCond),p.layer,p.nStimSets);
    
    fprintf('\n\nSESSION %d, RAT %d\n', sess, p.ratNum);
    
    %% execute model
    p = model(p,stims,weights);
    
    
    %% save it up
    p.runningTime = GetSecs-startTime;
    fName = strcat(p.dataDir, '/Session',num2str(sess),'_Rat',num2str(p.ratNum),'.mat');
    save(fName,'p');
    
end

fprintf ('\nFinished running rat %d \r', p.ratNum);
end