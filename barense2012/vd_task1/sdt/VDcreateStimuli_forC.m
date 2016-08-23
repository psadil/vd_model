function [p, stims] = VDcreateStimuli_forC(p)
% VDcreateStimuli - creates stimuli for input to VD network

% September 11, 2015 ps
% -> this VDcreateStimuli (unlike what gets created for when there are 5
% caudal grids) creates two dimensional inputs. Each dimension can take one
% of four different simple visual features (.1, .3, .7, .9). That makes 16
% possible 2-d inputs (4^2), and 4^8 (65536) 8-d inputs for the PRC layer

global consts

ROOT = [pwd '\'];
exptDir = [consts.exptName, '\'];

nTotalDims = p.components;
nCaudalGrids = p.numGrids_Caudal;
nDimsCaudal = nTotalDims / nCaudalGrids;
nStimFactors = p.nStimFactors;
nSimpleConj = nStimFactors ^ nDimsCaudal;



% number of features to use for stimuli (less then total)
nStimsPerCond = nSimpleConj;

LA_misMatch = zeros(p.nMismatch,nTotalDims,2); %one col for stim1, one col for stim2
LA_match = zeros(p.nMatch,nTotalDims,2);
HA_misMatch = zeros(p.nMismatch,nTotalDims,2); %one col for stim1, one col for stim2
HA_match = zeros(p.nMatch,nTotalDims,2);

% count=1;
% availFeat = zeros(nSimpleConj, nDimsCaudal);
% for inp1 = 1:nCaudalGrids,
%     for inp2 = 1:nStimFactors,
%         availFeat(count,:) = [inp1 inp2];
%         count=count+1;
%     end
% end
availFeat = permn(p.features,nDimsCaudal);

% availFeat(availFeat==1)=0;
% availFeat(availFeat==2)=1;

features = availFeat; %[repmat(availFeat(:,1),[1 3]) repmat(availFeat(:,2),[1 3])]; throwback to when each input dimension was replicated 3 times

firstFeatureToCheck=(1:nDimsCaudal:nTotalDims);
lastFeatureToCheck = firstFeatureToCheck+(nDimsCaudal-1);

notUnique = 1;

while notUnique
    
    %% LA
    
    
    %first create Mismatch pairs
    
    for stimPair = 1:p.nMismatch,
        %         chosenFeat=zeros(nTotalDims,2); %nDimsCaudal, 2 stimuli
        %         for feat = 1:nTotalDims, %'feature' in the sense of separate posterior grid
        %             pickFeatVal = randperm(8);
        %             chosenFeat(feat,:) = pickFeatVal(1:2); %this gives a different feature val for this feature across the 2 stimuli
        %         end
        chosenFeat = zeros(nCaudalGrids,2);
        while any(chosenFeat(:,1) == chosenFeat(:,2))
            chosenFeat = randi(nStimsPerCond,[nCaudalGrids,2]);
        end
        
        %         stimuli1(stimPair,:,1) = features(chosenFeat(:,1),:);
        %         stimuli2(stimPair,:,2) = features(chosenFeat(:,2),:);
        
        for feature = 1:nCaudalGrids
            %first stim in pair
            LA_misMatch(stimPair,firstFeatureToCheck(feature):lastFeatureToCheck(feature),1) = features(chosenFeat(feature,1),:);
            %second stim in pair
            LA_misMatch(stimPair,firstFeatureToCheck(feature):lastFeatureToCheck(feature),2) = features(chosenFeat(feature,2),:);
        end
        
    end
    
    %then create Match pairs
    for stimPair = 1:p.nMatch,
        
        chosenFeat = randi(nStimsPerCond,[nCaudalGrids,1]);
        for feature = 1:nCaudalGrids
            %first stim in pair
            LA_match(stimPair,firstFeatureToCheck(feature):lastFeatureToCheck(feature),1) = features(chosenFeat(feature,1),:);
        end
        
    end
    %second stim in pair (matches first)
    LA_match(:,:,2) = LA_match(:,:,1);
    
    stimuli1 = LA_misMatch;
    stimuli2 = LA_match;
    
    
    notUnique=0;
    breakOut = 0;
    for item = 1:p.nMismatch-1
        for trial = item+1:p.nMismatch
            if stimuli1(item,:,:) == (stimuli1(trial,:,:))
                notUnique = 1;
                breakOut = 1;
                break;
            elseif stimuli2(item,:,:) == (stimuli2(trial,:,:))
                notUnique=1;
                breakOut = 1;
                break;
            elseif stimuli1(item,:,:) == (stimuli2(trial,:,:))
                notUnique = 1;
                breakOut = 1;
                break;
            end
        end
        if breakOut
            break;
        end
    end
    
    %save stimuli
    if ~exist(strcat(ROOT,exptDir,'condition1'),'dir')
        mkdir(strcat(ROOT,exptDir,'condition1'))
    end
    location = strcat(ROOT,exptDir,'condition1', '/stimuli.mat');
    save(location,'stimuli1','stimuli2');
    
    
    %% high ambiguity
    
    for stimPair = 1:p.nMismatch,
        nEqual = 0;
        while nEqual ~= (nCaudalGrids - 1)
            nEqual = 0;
            chosenFeat = randi(nStimsPerCond,[nCaudalGrids,2]);
            
            for featInput = 1:nCaudalGrids
                stim1Feat = chosenFeat(featInput,1);
                stim2Feat = chosenFeat(featInput,2);
                if all(stim1Feat == stim2Feat)
                    nEqual = nEqual + 1;
                end
            end
            
            if nEqual == (nCaudalGrids - 1)
                %                 break
            else
                nEqual = 0;
            end
        end
        
        for grid = 1:nCaudalGrids
            %first stim in pair
            HA_misMatch(stimPair,firstFeatureToCheck(grid):lastFeatureToCheck(grid),1) = features(chosenFeat(grid,1),:);
            %second stim in pair
            HA_misMatch(stimPair,firstFeatureToCheck(grid):lastFeatureToCheck(grid),2) = features(chosenFeat(grid,2),:);
        end
        
        
    end
    
    
    %then create Match pairs
    for stimPair = 1:p.nMatch,
        
        chosenFeat = randi(nStimsPerCond,[nCaudalGrids,1]);
        for feature = 1:nCaudalGrids
            %first stim in pair
            HA_match(stimPair,firstFeatureToCheck(feature):lastFeatureToCheck(feature),1) = features(chosenFeat(feature,1),:);
        end
        
    end
    %second stim in pair (matches first)
    HA_match(:,:,2) = HA_match(:,:,1);
    
    stims.LA_misMatch = LA_misMatch;
    stims.LA_match = LA_match;
    stims.HA_misMatch = HA_misMatch;
    stims.HA_match = HA_match;
    
    stimuli1 = HA_misMatch;
    stimuli2 = HA_match;
    
    
    for item = 1:p.nMismatch-1
        
        for trial = item+1:p.nMismatch
            if stimuli1(item,:,:) == (stimuli1(trial,:,:))
                notUnique=1;
                breakOut = 1;
                break;
            elseif stimuli2(item,:,:) == (stimuli2(trial,:,:))
                notUnique=1;
                breakOut = 1;
                break;
            elseif stimuli1(item,:,:) == (stimuli2(trial,:,:))
                notUnique=1;
                breakOut = 1;
                break;
            end
        end
        if breakOut
            break;
        end
    end
    
    
    %save stimuli
    if ~exist(strcat(ROOT,exptDir,'condition2'),'dir')
        mkdir(strcat(ROOT,exptDir,'condition2'))
    end
    location = strcat(ROOT,exptDir,'condition2', '/stimuli.mat');
    save(location,'stimuli1','stimuli2');
    
end



