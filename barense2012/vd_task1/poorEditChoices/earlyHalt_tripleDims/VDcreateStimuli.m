function VDcreateStimuli(ambigCond)
% VDcreateStimuli - creates stimuli for input to VD network

% September 11, 2015 ps
% -> this VDcreateStimuli (unlike what gets created for when there are 5
% caudal grids) creates two dimensional inputs. Each dimension can take one
% of four different simple visual features (.1, .3, .7, .9). That makes 16
% possible 2-d inputs (4^2), and 4^8 (65536) 8-d inputs for the PRC layer

ROOT = [pwd '\'];
exptDir = 'september11_2015/';
nMatch = 36;
nMismatch= 36;
nTrials = nMismatch+nMatch;

nTotalDims = 8;
nSimpleFeats = 4;
nCaudalGrids = 4; % just chance that this equals nSimpleFeats here
nDimsCaudal = nTotalDims / nSimpleFeats;
nSimpleConj = nDimsCaudal ^ nSimpleFeats;


stimuli1 = zeros(nMismatch,nTotalDims,2); %one col for stim1, one col for stim2
stimuli2 = zeros(nMatch,nTotalDims,2);

count=1;
availFeat = zeros(nTotalDims, nDimsCaudal);
for inp1 = 1:nSimpleFeats,
    for inp2 = 1:nSimpleFeats,
        availFeat(count,:) = [inp1 inp2];
        count=count+1;
    end
end

availFeat(availFeat==1)=0.1;
availFeat(availFeat==2)=0.3;
availFeat(availFeat==3)=0.7;
availFeat(availFeat==4)=0.9;

features = availFeat; %[repmat(availFeat(:,1),[1 3]) repmat(availFeat(:,2),[1 3])]; throwback to when each input dimension was replicated 3 times

firstFeatureToCheck=(1:nDimsCaudal:nTotalDims);
lastFeatureToCheck = firstFeatureToCheck+(nDimsCaudal-1);

%first create Mismatch pairs
if ambigCond==1, %% low ambiguity, no features shared
    for stimPair = 1:nMismatch,
        %         chosenFeat=zeros(nTotalDims,2); %nDimsCaudal, 2 stimuli
        %         for feat = 1:nTotalDims, %'feature' in the sense of separate posterior grid
        %             pickFeatVal = randperm(8);
        %             chosenFeat(feat,:) = pickFeatVal(1:2); %this gives a different feature val for this feature across the 2 stimuli
        %         end
        chosenFeat = zeros(nCaudalGrids,2);
        while any(chosenFeat(:,1) == chosenFeat(:,2))
            chosenFeat = randi(nSimpleConj,[nCaudalGrids,2]);
        end
        
%         stimuli1(stimPair,:,1) = features(chosenFeat(:,1),:);
%         stimuli2(stimPair,:,2) = features(chosenFeat(:,2),:);
        
        for feature = 1:nCaudalGrids
            %first stim in pair
            stimuli1(stimPair,firstFeatureToCheck(feature):lastFeatureToCheck(feature),1) = features(chosenFeat(feature,1),:);
            %second stim in pair
            stimuli1(stimPair,firstFeatureToCheck(feature):lastFeatureToCheck(feature),2) = features(chosenFeat(feature,2),:);
        end
        
    end
    
    % high ambiguity
elseif ambigCond==2,
    for stimPair = 1:nMismatch,
        
        nEqual = 0;
        while nEqual ~= (nCaudalGrids - 1)
            nEqual = 0;
            chosenFeat = randi(nSimpleConj,[nCaudalGrids,2]);
            
            for featInput = 1:nCaudalGrids
                stim1Feat = chosenFeat(featInput,1);
                stim2Feat = chosenFeat(featInput,2);
                if all(stim1Feat == stim2Feat)
                    nEqual = nEqual + 1;
                end
            end
%             if nEqual > 1
%                 nEqual
%             end            
            if nEqual == (nCaudalGrids - 1) 
%                 break
                2;  
            else
                nEqual = 0;
            end
        end
                
%         stimuli1(stimPair,:,1) = features(chosenFeat(:,1),:);
%         stimuli2(stimPair,:,2) = features(chosenFeat(:,2),:);

        for grid = 1:nCaudalGrids
            %first stim in pair
            stimuli1(stimPair,firstFeatureToCheck(grid):lastFeatureToCheck(grid),1) = features(chosenFeat(grid,1),:);
            %second stim in pair
            stimuli1(stimPair,firstFeatureToCheck(grid):lastFeatureToCheck(grid),2) = features(chosenFeat(grid,2),:);
        end
        
        
        %
        %         chosenFeat=zeros(nDimsCaudal,2); % nDimsCaudal, 2 stimuli
        %         chooseMatchingFeat = randperm(3);
        %         chooseMatchingFeat = chooseMatchingFeat(1:2); % choose 2 features to match
        %         % Feat 1 is fill pattern, Feat 2 is inner shape, Feat 3 is outer shape
        %         % Feat 1 maps to feat1, Feat 2 maps to feats 2 and 3, Feat 3 maps to feats 4 and 5
        %         inpFeatMatch = zeros(3,1);
        %         inpFeatMatch(chooseMatchingFeat)=1;
        %         inpFeatMatch = cat(2,inpFeatMatch(1), inpFeatMatch(2),inpFeatMatch(2),inpFeatMatch(3),inpFeatMatch(3));
        %         for feat = 1:nDimsCaudal, %make all nDimsCaudal features
        %             pickFeatVal = randperm(nTotalDims);
        %             if inpFeatMatch(feat) == 1,
        %                 chosenFeat(feat,:) = repmat(pickFeatVal(1),[1 2]);
        %             elseif inpFeatMatch(feat) == 0,
        %                 chosenFeat(feat,:) = pickFeatVal(1:2);
        %             end
        %         end
        %
        %         %first stim in pair
        %         stimuli1(stimPair,1:3,1) = features(chosenFeat(1,1),:);
        %         stimuli1(stimPair,4:6,1) = features(chosenFeat(2,1),:);
        %         stimuli1(stimPair,7:9,1) = features(chosenFeat(3,1),:);
        %         stimuli1(stimPair,10:12,1) = features(chosenFeat(4,1),:);
        %         stimuli1(stimPair,13:15,1) = features(chosenFeat(5,1),:);
        %         %second stim in pair
        %         stimuli1(stimPair,1:3,2) = features(chosenFeat(1,2),:);
        %         stimuli1(stimPair,4:6,2) = features(chosenFeat(2,2),:);
        %         stimuli1(stimPair,7:9,2) = features(chosenFeat(3,2),:);
        %         stimuli1(stimPair,10:12,2) = features(chosenFeat(4,2),:);
        %         stimuli1(stimPair,13:15,2) = features(chosenFeat(5,2),:);
        
    end
end

        % number of stimuli that will have too few features identical
%         nWrongStims = sum((sum(stimuli1(:,:,1)==stimuli1(:,:,2),2)>=6)==0)


%then create Match pairs
for stimPair = 1:nMatch,
    
    chosenFeat = randi(nSimpleConj,[nCaudalGrids,1]);
    for feature = 1:nCaudalGrids
        %first stim in pair
        stimuli2(stimPair,firstFeatureToCheck(feature):lastFeatureToCheck(feature),1) = features(chosenFeat(feature,1),:);
    end
    
end
%second stim in pair (matches first)
stimuli2(:,:,2) = stimuli2(:,:,1);


%save stimuli
if ~exist(strcat(ROOT,exptDir,'condition', num2str(ambigCond)),'dir')
    mkdir(strcat(ROOT,exptDir,'condition', num2str(ambigCond)))
end
location = strcat(ROOT,exptDir,'condition', num2str(ambigCond), '/stimuli.mat');
save(location,'stimuli1','stimuli2');
