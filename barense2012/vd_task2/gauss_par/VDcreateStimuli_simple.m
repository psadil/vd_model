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
nStimsPerCond = {[10,16],[1,6]}; %[LA, HA]

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


% LA trials
cond = 1;
first1 = permn(nStimsPerCond{cond}(1):nStimsPerCond{cond}(2),p.numInputDims_PRC/p.numInputDims_Caudal);
final1 = zeros(size(first1,1),p.numInputDims_PRC);
firstIdx = 1;
for col = 1:p.numInputDims_Caudal:p.numInputDims_PRC
    final1(:,col) = first1(:,firstIdx);
    firstIdx=firstIdx+1;
end

for feature = nStimsPerCond{cond}(1):nStimsPerCond{cond}(2)
    idx = find(final1==feature);
    
    final1(idx)=availFeat(feature,1);
    final1(idx+size(first1,1))=availFeat(feature,2);    
end

% HA trials
cond = 2;
first2 = permn(nStimsPerCond{cond}(1):nStimsPerCond{cond}(2),p.numInputDims_PRC/p.numInputDims_Caudal);
final2 = zeros(size(first2,1),p.numInputDims_PRC);
firstIdx = 1;
for col = 1:p.numInputDims_Caudal:p.numInputDims_PRC
    final2(:,col) = first2(:,firstIdx);
    firstIdx=firstIdx+1;
end

for feature = nStimsPerCond{cond}(1):nStimsPerCond{cond}(2)
    idx = find(final2==feature);
    
    final2(idx)=availFeat(feature,1);
    final2(idx+size(first2,1))=availFeat(feature,2);    
end


nMatchingLA = 0;
tryAgain = 1;
while tryAgain
    breakOut = 0;
    
    % use these values to pull unique stimuli
    trials_LA = reshape(randsample(size(final1,1),p.nMismatch*2),[p.nMismatch,2]);
    trials_LA_misMatch1 = trials_LA(:,1);
    trials_LA_misMatch2 = zeros(size(trials_LA_misMatch1));
    
    % grab match for LA
    used = trials_LA(:); % the stims that we've already used
    for stim = 1:p.nMismatch
        if breakOut
            break
        end
        probeWith = first1(trials_LA_misMatch1(stim),:);
        
        row_randIdx = randperm(size(first1,1));
        for row = 1:size(first1,1)
            
            % don't try this row if we've already used it somewhere
            if ismember(row_randIdx(row),used)
                continue
            end
            
            probed = first1(row_randIdx(row),:);
            if sum(probeWith == probed) == nMatchingLA
                
                % grab this stim from the pre-tmp
                trials_LA_misMatch2(stim) = row_randIdx(row);
                
                % declare that row as used
                used = [used;row_randIdx(row)];
                
                % stop searching for a matching stim
                break % the row loop
            end
            
            % if we've reached the end of all rows
            if row == size(first1,1)
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
    trials_HA = reshape(randsample(size(final2,1),p.nMismatch*2),[p.nMismatch,2]);
    trials_HA_misMatch1 = trials_HA(:,1);
    trials_HA_misMatch2 = zeros(size(trials_HA_misMatch1));
    
    % grab match for LA
    used = trials_HA(:); % the stims that we've already used
    for stim = 1:p.nMismatch
        if breakOut
            break
        end
        probeWith = first2(trials_HA_misMatch1(stim),:);
        
        row_randIdx = randperm(size(first2,1));
        for row = 1:size(first2,1)
            
            % don't try this row if we've already used it somewhere
            if ismember(row_randIdx(row),used)
                continue
            end
            
            probed = first2(row_randIdx(row),:);
            if sum(probeWith == probed) == nMatchingHA
                
                % grab this stim from the pre-tmp
                trials_HA_misMatch2(stim) = row_randIdx(row);
                
                % declare that row as used
                used = [used;row_randIdx(row)];
                
                % stop searching for a matching stim
                break % the row loop
            end
            
            % if we've reached the end of all rows
            if row == size(first2,1)
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
stims.LA_match = zeros(p.nMismatch,8,2);
stims.LA_match(:,:,1) = final1(trials_LA(:,2),:);
stims.LA_match(:,:,2) = final1(trials_LA(:,2),:);

stims.HA_match = zeros(p.nMismatch,8,2);
stims.HA_match(:,:,1) = final2(trials_HA(:,2),:);
stims.HA_match(:,:,2) = final2(trials_HA(:,2),:);

stims.LA_misMatch = zeros(p.nMismatch,8,2);
stims.LA_misMatch(:,:,1) = final1(trials_LA_misMatch1,:);
stims.LA_misMatch(:,:,2) = final1(trials_LA_misMatch2,:);

stims.HA_misMatch = zeros(p.nMismatch,8,2);
stims.HA_misMatch(:,:,1) = final2(trials_HA_misMatch1,:);
stims.HA_misMatch(:,:,2) = final2(trials_HA_misMatch2,:);


end