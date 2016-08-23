function stims = init_stimuli(p)
% createListLengthStimuli - creates stimuli used for an experiment with a
% given rat. These stims are seen in every session. Stims are composed in
% pairs, such that one stim is presented for study while the other remains
% novel. Relative recognition of these two stimuli will eventually be
% compared.

% called by: runSim
% calls: permn

% input:
%   p: experimental structure. Used for many dimension and structure params

% output:
%   stims: array containing the stims that will be input.
%      nTrials x components x 2 x numListLengths x numStimSets.
%      NOTE: 2 refers to two stimuli per pair, one to be trained upon, one
%      that remains novel.

% NOTE: although every item is trial unique within a given
% stimSet/stimCond, stimuli are allowed to repeat across sets and
% conditions.


stims.misMatch = zeros(max(p.nTrials),p.components,2,length(p.nMismatch),p.nStimSets);
   
% make as many stims as stimSets requested
for stimSet = 1:p.nStimSets;
    
    % each stim condition has a different number of stims
    for stimCond = 1:length(p.nMismatch)
        
        % all of the stim generation happens in the following line. The
        % permn command generates everypossible permutation from which to
        % pull, and the encapsulating datasample command pulls (without
        % replacement) the number of needed stimuli.
        stimToPut = datasample(permn(p.features,p.components),p.nMismatch(stimCond)*2, 'Replace',false);
        
        % grab all of the trial unique stims for the trial-unique pairs
        stims.misMatch(1:p.nMismatch(stimCond),:,1,stimCond,stimSet) = ...
            stimToPut(1:(p.nMismatch(stimCond)),:);
        stims.misMatch(1:p.nMismatch(stimCond),:,2,stimCond,stimSet) = ...
            stimToPut(1+(p.nMismatch(stimCond)):end,:);
        
        % used in expt 3, when there are some matching trials
        if p.nMatch > 0
            
            stimToPut = datasample(permn([0,1],p.components),p.nMismatch(stimCond), 'Replace',false);

            % grab only the first and second stim generated from the repeating stims
            stims.match = zeros(size(stims.misMatch));
            stims.match(1:p.nMatch(stimCond),:,1,stimCond,stimSet) = stimToPut;
            stims.match(1:p.nMatch(stimCond),:,2,stimCond,stimSet) = stimToPut;
        end
        
        
    end
    
end

end