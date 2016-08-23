function [p, stims] = VDcreateStimuli_simple(p)
% VDcreateStimuli - creates stimuli for input to VD network

% September 11, 2015 ps
% -> this VDcreateStimuli (unlike what gets created for when there are 5
% caudal grids) creates two dimensional inputs. Each dimension can take one
% of four different simple visual features (.1, .3, .7, .9). That makes 16
% possible 2-d inputs (4^2), and 4^8 (65536) 8-d inputs for the PRC layer

% USE THIS FUNCTION TO CREATE STIMS!!!!!!!!!!!

nTotalDims = p.components;
nCaudalGrids = p.numGrids_Caudal;
nDimsCaudal = nTotalDims / nCaudalGrids;
nStimFactors = p.nStimFactors;
nSimpleConj = nDimsCaudal ^ nStimFactors;


% number of features to use for stimuli (less then total)
nStimsPerCond = 6;

% LA_misMatch = zeros(p.nMismatch,nTotalDims,2); %one col for stim1, one col for stim2
% LA_match = zeros(p.nMatch,nTotalDims,2);
% HA_misMatch = zeros(p.nMismatch,nTotalDims,2); %one col for stim1, one col for stim2
% HA_match = zeros(p.nMatch,nTotalDims,2);

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

nMatchingLA = 0;
tryAgain = 1;
while tryAgain
    breakOut = 0;
    
    % use these values to pull unique stimuli
    trials_LA = reshape(randsample(size(final,1),36*2),[36,2]);
    trials_LA_misMatch1 = trials_LA(:,1);
    trials_LA_misMatch2 = zeros(size(trials_LA_misMatch1));
    
    % grab match for LA
    used = trials_LA(:); % the stims that we've already used
    for stim = 1:36
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



%% now, make HA stims
nMatchingHA = 3;

tryAgain = 1;
while tryAgain
    breakOut = 0;
    
    % use these values to pull unique stimuli
    trials_HA = reshape(randsample(size(final,1),36*2),[36,2]);
    trials_HA_misMatch1 = trials_HA(:,1);
    trials_HA_misMatch2 = zeros(size(trials_HA_misMatch1));
    
    % grab match for LA
    used = trials_HA(:); % the stims that we've already used
    for stim = 1:36
        if breakOut
            break
        end
        probeWith = first(trials_HA_misMatch1(stim),:);
        
        row_randIdx = randperm(size(first,1));
        for row = 1:size(first,1)
            
            % don't try this row if we've already used it somewhere
            if ismember(row_randIdx(row),used)
                continue
            end
            
            probed = first(row_randIdx(row),:);
            if sum(probeWith == probed) == nMatchingHA
                
                % grab this stim from the pre-tmp
                trials_HA_misMatch2(stim) = row_randIdx(row);
                
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
stims.LA_match = zeros(36,8,2);
stims.LA_match(:,:,1) = final(trials_LA(:,2),:);
stims.LA_match(:,:,2) = final(trials_LA(:,2),:);

stims.HA_match = zeros(36,8,2);
stims.HA_match(:,:,1) = final(trials_HA(:,2),:);
stims.HA_match(:,:,2) = final(trials_HA(:,2),:);

stims.LA_misMatch = zeros(36,8,2);
stims.LA_misMatch(:,:,1) = final(trials_LA_misMatch1,:);
stims.LA_misMatch(:,:,2) = final(trials_LA_misMatch2,:);

stims.HA_misMatch = zeros(36,8,2);
stims.HA_misMatch(:,:,1) = final(trials_HA_misMatch1,:);
stims.HA_misMatch(:,:,2) = final(trials_HA_misMatch2,:);


% gather all stims for output, repeat as necessary
% stims.stimuli1 = reshape(repmat(reshape(LA_misMatch,...
%     [1,p.nMismatch,1,p.numInputDims_PRC/p.nDimReps,1,2]),[1,1,p.nDimReps,1,1,1]),...
%     [p.nMismatch,(p.numInputDims_PRC/p.nDimReps)*p.nDimReps,2]);

end