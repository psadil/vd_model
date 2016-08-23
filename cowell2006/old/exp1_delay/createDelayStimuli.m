function [p, stims] = createDelayStimuli(p)
% createDelayStimuli - creates stimuli for input to simulation of delay

% generates 'mis-matching' stims

nTotalDims = p.components;
nCaudalGrids = p.numGrids_Caudal;
nDimsCaudal = nTotalDims / nCaudalGrids;
nStimFactors = p.nStimFactors;
nSimpleConj = nDimsCaudal ^ nStimFactors;


% number of features to use for stimuli (less then total)
nStimsPerCond = 16;


count=1;
availFeat = zeros(nSimpleConj, nDimsCaudal);
for inp1 = 1:nCaudalGrids,
    for inp2 = 1:nStimFactors,
        availFeat(count,:) = [inp1 inp2];
        %     availFeat(count) = inp2;
        count=count+1;
    end
end

availFeat(availFeat==1)=0.05;
availFeat(availFeat==2)=0.35;
availFeat(availFeat==3)=0.65;
availFeat(availFeat==4)=0.95;

availFeat = Shuffle(availFeat,2);

first = permn(1:nStimsPerCond,p.numInputDims_PRC/p.numInputDims_Caudal);
final = zeros(size(first,1),p.numInputDims_PRC);
firstIdx = 1;
for col = 1:p.numInputDims_Caudal:p.numInputDims_PRC
    final(:,col) = first(:,firstIdx);
    firstIdx=firstIdx+1;
end

for feature = 1:nStimsPerCond
    idx = find(final==feature);
    
    final(idx)=availFeat(feature,1);
    final(idx+size(first,1))=availFeat(feature,2);
    
end

%%
stims = zeros(max(p.nTrials), p.components,2);

nMatchingLA = 0;
tryAgain = 1;
while tryAgain
    breakOut = 0;
    
    % use these values to pull unique stimuli
    trials_LA = reshape(randsample(size(final,1),p.nMismatch),...
        [p.nMismatch,1]);
    trials_LA_misMatch1 = trials_LA(:,1);
    trials_LA_misMatch2 = zeros(size(trials_LA_misMatch1));
    
    % grab match for LA
    used = trials_LA(:); % the stims that we've already used
    for stim = 1:p.nMismatch
        if breakOut
            break
        end
        probeWith = first(trials_LA_misMatch1(stim),:);
        
        row_randIdx = randperm(size(first,1));
        for row = 1:size(first,1)
            
            % don't try this row if we've already used it somewhere
            if ismember(row_randIdx(row),used)
                continue
            end
            
            probed = first(row_randIdx(row),:);
            if sum(probeWith == probed) == nMatchingLA
                
                % grab this stim from the pre-tmp
                trials_LA_misMatch2(stim) = row_randIdx(row);
                
                % declare that row as used
                used = [used;row_randIdx(row)];
                
                % stop searching for a matching stim
                break % the row loop
            end
            
            % if we've reached the end of all rows
            if row == size(first,1)
                breakOut = 1;
            end
        end
        
        % for when we broke out of the rows
        if breakOut
            break
        end
        
    end
    
    % for when we've broken out of the stim loop
    if breakOut
        continue
    end
    
    % but, if we've made it this far..
    tryAgain=0;
    
end


%% classify as features

% grab all of the trial unique stims for the trial-unique pairs
stims(:,:,1) = final(trials_LA_misMatch1,:);
stims(:,:,2) = final(trials_LA_misMatch2,:);



end