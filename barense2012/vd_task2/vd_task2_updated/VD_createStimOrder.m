function [ p ] = VD_createStimOrder( p )
%VD_createStimOder creates stimulus order to be presented to simulation
%during a given session
% called during run_sim

% October 31, 2015 ps
% NOW, creates stimuli for a different task, 88 trials long.

% HA block -> 88 HA objects (44 match, 44 non-match)
% LA block -> 30 HA objects (15 match, 15 non-match), 58 LA objects (29 match, 29 misMatch)

% HA object refers to stimuli pulled from 1-6 (and share 4 features on
% non-match trials, LA objects are pulled from features 7-16 (and share no
% features on misMatching trials)

trials = 1:p.nTrials;
compTrials = 1:3:p.nTrials;
trials(compTrials)=[];


%check to see whether there are more than 3 trials in a row the same % if yes, redo
check = 1;
while check
    p.tType = cat(2,ones(1,p.nMismatch), 2*ones(1,p.nMatch));
    p.tType = p.tType(randperm(length(p.tType)));
    
    
    % first, check to make sure that there are an evenly distributed amount of
    % match (2) and misMatch (1) trials on the compTrials
    nMatching = sum(p.tType(compTrials)==2);
    nMisMatching = sum(p.tType(compTrials)==1);
    
    
    if nMatching== p.nComps/2 && nMisMatching == p.nComps/2
        check = 0;
    end
    
    %% then -- if in an LA trial -- that the intervening trials (trials that aren't
    % in compTrials) are from LA stims, while the rest aren't
    
    if any(p.stimCond == [1,3])
        p.tType(trials) = p.tType(trials)+2;
    end
    
    
    p.stimOrder = [randperm(p.nMismatch)' randperm(p.nMatch)'];
    
end

