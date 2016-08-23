function [] = plotFamilDiffs(firstRat, lastRat, folderName, onlyFigure)
% analyze familiarity differences

% now called directly from create_sim

% currently, taking the maximum familiarity differences at each trial
% produces the set of most desirable results
% folderName must be a string (without / on either side)

% folderName = '1encod_p6Ratio1p2_1sampVar5_0noise_100train_20Max25_5peak_100rows_A01_B4_forceSamp_NOshrinkingLearn_etaExp5p_Gexpp8_stim4Diff_NOfives';

% scrsz = get(groot, 'ScreenSize');

nSess = 6;
layers = 2;
nSessPerLayer = nSess/layers;

saveFolder = [pwd,'/graphsAndSession/', folderName];

% load sample data file to get size of parameters
fileName = [saveFolder, '/Session', num2str(1), '_Rat', num2str(1)];
load(fileName)

numRats = lastRat-firstRat+1;

tType = zeros(numRats,nSess,p.nTrials);

familDifferences_misMatch_caudal = zeros(numRats,nSess,p.nTrials);
familDifferences_match_caudal = zeros(numRats,nSess,p.nTrials);
familDifferences_match_PRC = zeros(numRats,nSess,p.nTrials);
familDifferences_misMatch_PRC = zeros(numRats,nSess,p.nTrials);

% familDiff_caudal = zeros(numRats,nSess,p.nTrials);
% familDiff_PRC = zeros(numRats,nSess,p.nTrials);

answer = zeros(numRats,nSess,p.nTrials);
correct = zeros(numRats,nSess,p.nTrials);
acc_firstHalf = zeros(numRats,nSess);
acc_secondHalf = zeros(numRats,nSess);
acc_match_first = zeros(numRats,nSess);
acc_misMatch_first = zeros(numRats,nSess);
acc_match_second = zeros(numRats,nSess);
acc_misMatch_second = zeros(numRats,nSess);

dPrime_first = zeros(numRats,nSess);
dPrime_second = zeros(numRats,nSess);

hitRate_first = zeros(numRats,nSess);
hitRate_second = zeros(numRats,nSess);
FARate_first = zeros(numRats,nSess);
FARate_second = zeros(numRats,nSess);

meanSelectivity_caudal_prev = zeros(numRats,nSess,p.nTrials);
meanSelectivity_caudal_new = zeros(numRats,nSess,p.nTrials);
meanSelectivity_PRC_prev = zeros(numRats,nSess,p.nTrials);
meanSelectivity_PRC_new = zeros(numRats,nSess,p.nTrials);

familDiff_used = zeros(numRats,nSess,p.nTrials);

meanSelectivity_caudal_prev_misMatch = zeros(numRats,nSess,p.nTrials);
meanSelectivity_caudal_prev_match = zeros(numRats,nSess,p.nTrials);
meanSelectivity_caudal_new_misMatch = zeros(numRats,nSess,p.nTrials);
meanSelectivity_caudal_new_match = zeros(numRats,nSess,p.nTrials);

meanSelectivity_PRC_prev_misMatch = zeros(numRats,nSess,p.nTrials);
meanSelectivity_PRC_prev_match = zeros(numRats,nSess,p.nTrials);
meanSelectivity_PRC_new_misMatch = zeros(numRats,nSess,p.nTrials);
meanSelectivity_PRC_new_match = zeros(numRats,nSess,p.nTrials);

comparedFeat = zeros(numRats,nSess,p.nTrials,p.numGrids_Caudal);

trial_misMatch_caudal = zeros(nSess,p.nTrials);
trial_match_caudal = zeros(nSess,p.nTrials);
trial_misMatch_PRC = zeros(nSess,p.nTrials);
trial_match_PRC = zeros(nSess,p.nTrials);

trial_misMatch_caudal_if = zeros(numRats,nSess,p.nTrials);
trial_match_caudal_if = zeros(numRats,nSess,p.nTrials);
trial_misMatch_PRC_if = zeros(numRats,nSess,p.nTrials);
trial_match_PRC_if = zeros(numRats,nSess,p.nTrials);

gridUsed_one_total = zeros(nSess,p.nTrials);
gridUsed_one = zeros(numRats,nSess,p.nTrials);

gridUsed_one_misMatch_total = zeros(nSess,p.nTrials);
gridUsed_one_misMatch = zeros(numRats,nSess,p.nTrials);

gridUsed_one_match_total = zeros(nSess,p.nTrials);
gridUsed_one_match = zeros(numRats,nSess,p.nTrials);

selectivity_caudal_prev_match_grid1 = zeros(numRats,nSess,p.nTrials);
selectivity_caudal_prev_misMatch_grid1 = zeros(numRats,nSess,p.nTrials);
selectivity_caudal_new_match_grid1 = zeros(numRats,nSess,p.nTrials);
selectivity_caudal_new_misMatch_grid1 = zeros(numRats,nSess,p.nTrials);


fixationByComparison_prev = zeros(numRats,nSess,p.nTrials);
fixationByComparison_new = zeros(numRats,nSess,p.nTrials);

fixByComp_prev_mean = zeros(numRats,nSess,p.nTrials,max(p.maxFixations));
fixByComp_new_mean = zeros(numRats,nSess,p.nTrials,max(p.maxFixations));


peakAct = zeros(numRats,nSess,p.nTrials,2);
totalAct = zeros(numRats,nSess,p.nTrials,2);

featSamplePerTrial = zeros(numRats,nSess,p.nTrials);
featSamplePerTrial_misMatch = zeros(numRats,nSess,p.nTrials);
featSamplePerTrial_match = zeros(numRats,nSess,p.nTrials);
featSamplePerTrial_misMatch_count = zeros(nSess,p.nTrials);
featSamplePerTrial_match_count = zeros(nSess,p.nTrials);

ifPRC_both = zeros(numRats,nSess,2,p.nTrials);

act_peak_prev_init = zeros(numRats,nSess,2,p.nTrials);
act_peak_prev_fin = zeros(numRats,nSess,2,p.nTrials);
act_peak_new_init = zeros(numRats,nSess,2,p.nTrials);
act_peak_new_fin = zeros(numRats,nSess,2,p.nTrials);
act_total_prev_init = zeros(numRats,nSess,2,p.nTrials);
act_total_prev_fin = zeros(numRats,nSess,2,p.nTrials);
act_total_new_init = zeros(numRats,nSess,2,p.nTrials);
act_total_new_fin = zeros(numRats,nSess,2,p.nTrials);

fixations_total = zeros(numRats,nSess,p.nTrials);

compTrials = 1:3:p.nTrials;

for rat = firstRat:lastRat
    for session = 1:nSess
        
        fileName = [saveFolder, '/Session', num2str(session), '_Rat', num2str(rat)];
        load(fileName)
        fprintf ('\nloading rat %d, session %d', rat, session);
        
        tType(rat,session,:) = p.tType;
        
        tType(tType==3) = 1;
        tType(tType==4) = 2;
        
        p.tType(p.tType==3) = 1;
        p.tType(p.tType==4) = 2;
        
        
        answer(rat,session,:) = p.answer;
        correct(rat,session,:) = p.correct;
        familDiff_used(rat,session,:) = p.famil_difference;
        comparedFeat(rat,session,:,:) = p.comparedFeat;
        fixations_total(rat,session,:) = p.fixations;
        
        %------------------------------------------------------------------
        % collect all trial-wise selectivity of choice phase
        %------------------------------------------------------------------
        
        % overall selectivity
        %         familDiff_caudal(rat,session,:) = p.familDiff_caudal;
        meanSelectivity_caudal_prev(rat,session,:) = p.meanSelectivity_caudal_prev;
        meanSelectivity_caudal_new(rat,session,:) = p.meanSelectivity_caudal_new;
        
        
        % if session will contain PRC layer, grab that layer too
        if p.which_gp_layer == 2
            %             familDiff_PRC(rat,session,:) = p.famil_difference;
            meanSelectivity_PRC_prev(rat,session,:) = p.meanSelectivity_PRC_prev;
            meanSelectivity_PRC_new(rat,session,:) = p.meanSelectivity_PRC_new;
        end
        
        %------------------------------------------------------------------
        % collect all trial-wise peak and total activations
        %------------------------------------------------------------------
        
        act_peak_prev_init(rat,session,1,:) = p.prevStimInit_act_peak(1,:);
        act_peak_prev_fin(rat,session,1,:) = p.prevStimFin_act_peak(1,:);
        act_peak_new_init(rat,session,1,:) = p.newStimInit_act_peak(1,:);
        act_peak_new_fin(rat,session,1,:) = p.newStimFin_act_peak(1,:);
        act_total_prev_init(rat,session,1,:) = p.prevStimInit_act_total(1,:);
        act_total_prev_fin(rat,session,1,:) = p.prevStimFin_act_total(1,:);
        act_total_new_init(rat,session,1,:) = p.newStimInit_act_total(1,:);
        act_total_new_fin(rat,session,1,:) = p.newStimFin_act_total(1,:);
        
        if p.which_gp_layer == 2
            
            act_peak_prev_init(rat,session,2,:) = p.prevStimInit_act_peak(2,:);
            act_peak_prev_fin(rat,session,2,:) = p.prevStimFin_act_peak(2,:);
            act_peak_new_init(rat,session,2,:) = p.newStimInit_act_peak(2,:);
            act_peak_new_fin(rat,session,2,:) = p.newStimFin_act_peak(2,:);
            act_total_prev_init(rat,session,2,:) = p.prevStimInit_act_total(2,:);
            act_total_prev_fin(rat,session,2,:) = p.prevStimFin_act_total(2,:);
            act_total_new_init(rat,session,2,:) = p.newStimInit_act_total(2,:);
            act_total_new_fin(rat,session,2,:) = p.newStimFin_act_total(2,:);
            
        end
        
        
        % need to also gather for these trials when PRC was actually used
        ifPRC_both(rat,session,:,:) = p.usePRC;
        
        
        if any(session == [1,3,4,6])
            whichFix = 1;
        else
            whichFix = 2;
        end
        
        
        %------------------------------------------------------------------
        % delve in to trial-wise stuff
        %------------------------------------------------------------------
        
        for trial = 1:p.nTrials
            
            
            %------------------------------------------------------------------
            % grab number of features sampled in each comparison
            %------------------------------------------------------------------
            
            % NOTE: this divides by stimulus number, NOT actually new/prev
            count_prev = 0;
            count_new = 0;
            temp_prev = 0;
            temp_new = 0;
            temp_per_count_prev = zeros(1,p.maxFixations(whichFix));
            count_per_prev = zeros(1,p.maxFixations(whichFix));
            temp_per_count_new = zeros(1,p.maxFixations(whichFix));
            count_per_new = zeros(1,p.maxFixations(whichFix));
            for comparison = 1:p.maxFixations(whichFix)
                
                if p.featsSampedByComparison(1,trial,comparison);
                    temp_prev = temp_prev + p.featsSampedByComparison(1,trial,comparison);
                    count_prev = count_prev + 1;
                    
                    temp_per_count_prev(comparison) = p.featsSampedByComparison(1,trial,comparison);
                    count_per_prev(comparison) = count_per_prev(comparison) + 1;
                end
                if p.featsSampedByComparison(2,trial,comparison);
                    temp_new = temp_new + p.featsSampedByComparison(2,trial,comparison);
                    count_new = count_new + 1;
                    
                    temp_per_count_new(comparison) = p.featsSampedByComparison(1,trial,comparison);
                    count_per_new(comparison) = count_per_prev(comparison) + 1;
                    
                end
                
            end
            
            %--------------------------------------------------------------
            % gather those counts and average
            %--------------------------------------------------------------
            
            fixationByComparison_prev(rat,session,trial) = temp_prev/count_prev;
            fixationByComparison_new(rat,session,trial) = temp_new/count_new;
            
            fixByComp_prev_mean(rat,session,trial,1:p.maxFixations(whichFix)) = temp_per_count_prev./count_per_prev;
            fixByComp_new_mean(rat,session,trial,1:p.maxFixations(whichFix)) = temp_per_count_new./count_per_new;
            
            featSamplePerTrial(rat,session,trial) = sum(temp_prev) + sum(temp_new);
            %--------------------------------------------------------------
            
            
            
            
            
            for grid = 1:p.numGrids_Caudal
                if p.comparedFeat(trial,grid)
                    %                     selecByGrid_caudal_prev(rat,session,grid) = p.meanSelectivity_caudal_prev(trial,1);
                    gridUsed_one(rat,session,trial) = 1;
                    gridUsed_one_total(session,trial) = gridUsed_one_total(session,trial) + 1;
                end
            end
            
            
            
            %--------------------------------------------------------------
            % misMatch Trials
            %--------------------------------------------------------------
            
            if any(p.tType(trial) == [1,3])
                
                %----------------------------------------------------------
                % grab number of features sampled in each comparison
                %----------------------------------------------------------
                
                % NOTE: this divides by stimulus number, NOT actually new/prev
                count_prev_misMatch = 0;
                count_new_misMatch = 0;
                temp_prev_misMatch = 0;
                temp_new_misMatch = 0;
                for comparison = 1:p.maxFixations(whichFix)
                    
                    if p.featsSampedByComparison(1,trial,comparison);
                        temp_prev_misMatch = temp_prev_misMatch + p.featsSampedByComparison(1,trial,comparison);
                        count_prev_misMatch = count_prev_misMatch + 1;
                    end
                    if p.featsSampedByComparison(2,trial,comparison);
                        temp_new_misMatch = temp_new_misMatch + p.featsSampedByComparison(2,trial,comparison);
                        count_new_misMatch = count_new_misMatch + 1;
                    end
                end
                featSamplePerTrial_misMatch(rat,session,trial) = sum(temp_prev_misMatch) + sum(temp_new_misMatch);
                featSamplePerTrial_misMatch_count(session,trial) = featSamplePerTrial_misMatch_count(session,trial) + 1;
                
                %----------------------------------------------------------
                % Gather selectivity for misMatching Trials
                %----------------------------------------------------------
                
                % caudal layer
                %                 meanSelectivity_caudal_prev_misMatch(rat,session,trial,:) = p.meanSelectivity_caudal_prev(trial);
                %                 meanSelectivity_caudal_new_misMatch(rat,session,trial,:) = p.meanSelectivity_caudal_new(trial);
                
                familDifferences_misMatch_caudal(rat,session,trial) = p.familDiff_caudal(trial);
                trial_misMatch_caudal(session,trial) = trial_misMatch_caudal(session,trial) + 1;
                trial_misMatch_caudal_if(rat,session,trial) = 1;
                
                % look at selectivity
                
                
                
                
                if p.which_gp_layer == 2 && (all(p.usePRC(:,trial)))
                    %                     meanSelectivity_PRC_prev_misMatch(rat,session,trial) = p.meanSelectivity_PRC_prev(trial);
                    %                     meanSelectivity_PRC_new_misMatch(rat,session,trial) = p.meanSelectivity_PRC_new(trial);
                    
                    familDifferences_misMatch_PRC(rat,session,trial) = p.familDiff_PRC(trial);
                    trial_misMatch_PRC(session,trial) = trial_misMatch_PRC(session,trial) + 1;
                    trial_misMatch_PRC_if(rat,session,trial) = 1;
                end
                
                %--------------------------------------------------------------
                % match Trials
                %--------------------------------------------------------------
                
            elseif any(p.tType(trial) == [2,4])
                
                
                % NOTE: this divides by stimulus number, NOT actually new/prev
                count_prev_match = 0;
                count_new_match = 0;
                temp_prev_match = 0;
                temp_new_match = 0;
                for comparison = 1:p.maxFixations(whichFix)
                    
                    if p.featsSampedByComparison(1,trial,comparison);
                        temp_prev_match = temp_prev_match + p.featsSampedByComparison(1,trial,comparison);
                        count_prev_match = count_prev_match + 1;
                    end
                    if p.featsSampedByComparison(2,trial,comparison);
                        temp_new_match = temp_new_match + p.featsSampedByComparison(2,trial,comparison);
                        count_new_match = count_new_match + 1;
                    end
                end
                
                %----------------------------------------------------------
                % gather those counts and average
                %----------------------------------------------------------
                featSamplePerTrial_match(rat,session,trial) = sum(temp_prev_match) + sum(temp_new_match);
                featSamplePerTrial_match_count(session,trial) = featSamplePerTrial_match_count(session,trial) + 1;
                
                
                %                 meanSelectivity_caudal_prev_match(rat,session,trial,:) = p.meanSelectivity_caudal_prev(trial,:);
                %                 meanSelectivity_caudal_new_match(rat,session,trial,:) = p.meanSelectivity_caudal_new(trial,:);
                
                familDifferences_match_caudal(rat,session,trial) = p.familDiff_caudal(trial);
                trial_match_caudal(session,trial) = trial_match_caudal(session,trial) + 1;
                trial_match_caudal_if(rat,session,trial) = 1;
                %----------------------------------------------------------
                
                
                
                % add PRC to the mix
                if p.which_gp_layer == 2 && (all(p.usePRC(:,trial)))
                    %                     meanSelectivity_PRC_prev_match(rat,session,trial) = p.meanSelectivity_PRC_prev(trial);
                    %                     meanSelectivity_PRC_new_match(rat,session,trial) = p.meanSelectivity_PRC_new(trial);
                    familDifferences_match_PRC(rat,session,trial) = p.familDiff_PRC(trial);
                    trial_match_PRC(session,trial) = trial_match_PRC(session,trial) + 1;
                    trial_match_PRC_if(rat,session,trial) = 1;
                    
                end
            end
            
        end
        
        %------------------------------------------------------------------
        % tally results
        %------------------------------------------------------------------
        
        hitRate_first(rat,session) = sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==1))/sum(p.tType(1:p.nTrials/2)==1);  % a 'yes' (mismatch judgement) on trials that were mismatches (ie, p.tType==1)
        FARate_first(rat,session) = sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==2))/sum(p.tType(1:p.nTrials/2)==2); % a 'yes' on matching trials
        
        hitRate_second(rat,session) = sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==1))/sum(p.tType(p.nTrials/2+1:end)==1);  % a 'yes' (mismatch judgement) on trials that were mismatches (ie, p.tType==1)
        FARate_second(rat,session) = sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==2))/sum(p.tType(p.nTrials/2+1:end)==2); % a 'yes' on matching trials
        
        % hit/FA rate, over the entire block
        hitRate(rat,session) = sum(p.answer.*(p.tType==1))/sum(p.tType==1);
        FARate(rat,session) = sum(p.answer.*(p.tType==2))/sum(p.tType==2);
        
        hitRate_comp(rat,session) = sum(p.answer(compTrials).*(p.tType(compTrials)==1))/sum(p.tType(compTrials)==1);
        FARate_comp(rat,session) = sum(p.answer(compTrials).*(p.tType(compTrials)==2))/sum(p.tType(compTrials)==2);
        
        % adjust by adding to all
        hitRate_first_adj_all(rat,session) = (sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==1)) + 1/(2*sum((p.tType(1:p.nTrials/2)==1))))/(sum(p.tType(1:p.nTrials/2)==1)+1);
        FARate_first_adj_all(rat,session) = (sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==2)) + 1/(2*sum((p.tType(1:p.nTrials/2)==2))))/(sum(p.tType(1:p.nTrials/2)==2)+1);
        
        hitRate_second_adj_all(rat,session) = (sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==1)) + 1/(2*sum((p.tType(1:p.nTrials/2)==1))))/(sum(p.tType(p.nTrials/2+1:end)==1)+1);  % a 'yes' (mismatch judgement) on trials that were mismatches (ie, p.tType==1)
        FARate_second_adj_all(rat,session) = (sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==2)) + 1/(2*sum((p.tType(1:p.nTrials/2)==2))))/(sum(p.tType(p.nTrials/2+1:end)==2)+1); % a 'yes' on matching trials
        
        % hit/FA rate, over the entire block
        hitRate_all(rat,session) = (sum(p.answer.*(p.tType==1)) + 1/(2*sum((p.tType==1))))/(sum(p.tType==1)+1);
        FARate_all(rat,session) = (sum(p.answer.*(p.tType==2)) + 1/(2*sum((p.tType==2))))/(sum(p.tType==2)+1);
        
        hitRate_all_comp(rat,session) = (sum(p.answer(compTrials).*(p.tType(compTrials)==1)) + 1/(2*sum((p.tType(compTrials)==1))))/(sum(p.tType(compTrials)==1)+1);
        FARate_all_comp(rat,session) = (sum(p.answer(compTrials).*(p.tType(compTrials)==2)) + 1/(2*sum((p.tType(compTrials)==2))))/(sum(p.tType(compTrials)==2)+1);
        
        
        % adjust by adding to only trials
        hitRate_first_adj_some(rat,session) = hitRate_first(rat,session);
        if hitRate_first_adj_some(rat,session) == 1
            hitRate_first_adj_some(rat,session) = hitRate_first_adj_some(rat,session) - 1/(2*sum((p.tType(1:p.nTrials/2)==1)));
        elseif hitRate_first_adj_some(rat,session) == 0
            hitRate_first_adj_some(rat,session) = hitRate_first_adj_some(rat,session) + 1/(2*sum((p.tType(1:p.nTrials/2)==1)));
        end
        
        FARate_first_adj_some(rat,session) = FARate_first(rat,session);
        if FARate_first_adj_some(rat,session) == 1
            FARate_first_adj_some(rat,session) = FARate_first_adj_some(rat,session) - 1/(2*sum((p.tType(1:p.nTrials/2)==2)));
        elseif FARate_first_adj_some(rat,session) == 0
            FARate_first_adj_some(rat,session) = FARate_first_adj_some(rat,session) + 1/(2*sum((p.tType(1:p.nTrials/2)==2)));
        end
        
        hitRate_second_adj_some(rat,session) = hitRate_second(rat,session);
        if hitRate_second_adj_some(rat,session) == 1
            hitRate_second_adj_some(rat,session) = hitRate_second_adj_some(rat,session) - 1/(2*sum((p.tType(1:p.nTrials/2)==1)));
        elseif hitRate_second_adj_some(rat,session) == 0
            hitRate_second_adj_some(rat,session) = hitRate_second_adj_some(rat,session) + 1/(2*sum((p.tType(1:p.nTrials/2)==1)));
        end
        
        FARate_second_adj_some(rat,session) = FARate_second(rat,session);
        if FARate_second_adj_some(rat,session) == 1
            FARate_second_adj_some(rat,session) = FARate_second_adj_some(rat,session) - 1/(2*sum((p.tType(1:p.nTrials/2)==2)));
        elseif FARate_second_adj_some(rat,session) == 0
            FARate_second_adj_some(rat,session) = FARate_second_adj_some(rat,session) + 1/(2*sum((p.tType(1:p.nTrials/2)==2)));
        end
        
        
        acc_firstHalf(rat,session) = p.Acc_firstHalf;
        acc_secondHalf(rat,session) = p.Acc_secondHalf;
        
        acc_match_first(rat,session) = sum(answer(rat,session,find(tType(rat,session,(1:p.nMismatch))==2))==0)/(sum(tType(rat,session,1:p.nMismatch)==2));
        acc_match_second(rat,session) = sum(answer(rat,session,find(tType(rat,session,(p.nMismatch+1:p.nTrials))==2))==0)/(sum(tType(rat,session,p.nMismatch+1:p.nTrials)==2));
        acc_misMatch_first(rat,session) = sum(answer(rat,session,find(tType(rat,session,(1:p.nMismatch))==1))==1)/(sum(tType(rat,session,1:p.nMismatch)==1));
        acc_misMatch_second(rat,session) = sum(answer(rat,session,find(tType(rat,session,(p.nMismatch+1:p.nTrials))==1))==1)/(sum(tType(rat,session,p.nMismatch+1:p.nTrials)==1));
        
        
        %         % divide selectivity by peak activation
        %         peakAct(rat,session,:,:) = p.peak_act;
        %         totalAct(rat,session,:,:) = p.totalAct;
    end
end

%%

fixByComp_prev = squeeze(mean(fixationByComparison_prev,1));
fixByComp_new = squeeze(mean(fixationByComparison_new,1));

fixByComp_prev_mean_plot = squeeze(mean(fixByComp_prev_mean,1));
fixByComp_new_mean_plot = squeeze(mean(fixByComp_new_mean,1));

featSamplePerTrial_mean = squeeze(mean(featSamplePerTrial,1));
%%
featSamplePerTrial_misMatch_temp = squeeze(sum(featSamplePerTrial_misMatch,1));
featSamplePerTrial_misMatch_mean = featSamplePerTrial_misMatch_temp ./ featSamplePerTrial_misMatch_count;

featSamplePerTrial_match_temp = squeeze(sum(featSamplePerTrial_match,1));
featSamplePerTrial_match_mean = featSamplePerTrial_match_temp ./ featSamplePerTrial_match_count;


%% peak activation, by tType
%(rat x session x layer x trial)

tTypeCount_misMatch = squeeze(sum(tType==1,1));
tTypeCount_match = squeeze(sum(tType==2,1));

ifPRC_used = squeeze(sum(ifPRC_both,3))==2;
ifCaudal_used = squeeze(sum(ifPRC_both,3))<2;


tTypeCount_misMatch_caudal = squeeze(sum(tType==1.*ifCaudal_used,1));
tTypeCount_match_caudal = squeeze(sum(tType==2.*ifCaudal_used,1));

tTypeCount_misMatch_PRC = squeeze(sum(tType==1.*ifPRC_used,1));
tTypeCount_match_PRC = squeeze(sum(tType==2.*ifPRC_used,1));

% *************************************************************************
% new
% *************************************************************************

%--------------------------------------------------------------------------
% peak
%--------------------------------------------------------------------------

% caudal (only looking at trials in which caudal was used as judge)
peakAct_misMatch_caudal_new_temp = squeeze(act_peak_new_init(:,:,1,:)).*(tType==1).*ifCaudal_used;
peakAct_misMatch_caudal_new_temp_count = squeeze(sum(peakAct_misMatch_caudal_new_temp,1));

peakAct_misMatch_caudal_new_avg = peakAct_misMatch_caudal_new_temp_count./tTypeCount_misMatch_caudal;

peakAct_match_caudal_new_temp = squeeze(act_peak_new_init(:,:,1,:)).*(tType==2).*ifCaudal_used;
peakAct_match_caudal_new_temp_count = squeeze(sum(peakAct_match_caudal_new_temp,1));

peakAct_match_caudal_new_avg = peakAct_match_caudal_new_temp_count./tTypeCount_match_caudal;

% PRC
peakAct_misMatch_PRC_temp_new = squeeze(act_peak_new_init(:,:,2,:)).*(tType==1).*ifPRC_used;
peakAct_misMatch_PRC_new_temp_count = squeeze(sum(peakAct_misMatch_PRC_temp_new,1));

peakAct_misMatch_PRC_new_avg = peakAct_misMatch_PRC_new_temp_count./tTypeCount_misMatch_PRC;

peakAct_match_PRC_new_temp = squeeze(act_peak_new_init(:,:,2,:)).*(tType==2).*ifPRC_used;
peakAct_match_PRC__new_temp_count = squeeze(sum(peakAct_match_PRC_new_temp,1));

peakAct_match_PRC_new_avg = peakAct_match_PRC__new_temp_count./tTypeCount_match_PRC;

%--------------------------------------------------------------------------
%total
%--------------------------------------------------------------------------

% caudal
totalAct_misMatch_caudal_new_temp = squeeze(act_total_new_init(:,:,1,:)).*(tType==1).*ifCaudal_used;
totalAct_misMatch_caudal_new_temp_count = squeeze(sum(totalAct_misMatch_caudal_new_temp,1));

totalAct_misMatch_caudal_new_avg = totalAct_misMatch_caudal_new_temp_count./tTypeCount_misMatch_caudal;

totalAct_match_caudal_new_temp = squeeze(act_total_new_init(:,:,1,:)).*(tType==2).*ifCaudal_used;
totalAct_match_caudal_new_temp_count = squeeze(sum(totalAct_match_caudal_new_temp,1));

totalAct_match_caudal_new_avg = totalAct_match_caudal_new_temp_count./tTypeCount_match_caudal;

% PRC
totalAct_misMatch_PRC_new_temp = squeeze(act_total_new_init(:,:,2,:)).*(tType==1).*ifPRC_used;
totalAct_misMatch_PRC_new_temp_count = squeeze(sum(totalAct_misMatch_PRC_new_temp,1));

totalAct_misMatch_PRC_new_avg = totalAct_misMatch_PRC_new_temp_count./tTypeCount_misMatch_PRC;

totalAct_match_PRC_new_temp = squeeze(act_total_new_init(:,:,2,:)).*(tType==2).*ifPRC_used;
totalAct_match_PRC_new_temp_count = squeeze(sum(totalAct_match_PRC_new_temp,1));

totalAct_match_PRC_new_avg = totalAct_match_PRC_new_temp_count./tTypeCount_match_PRC;

% *************************************************************************
% prev
% *************************************************************************

%--------------------------------------------------------------------------
% peak
%--------------------------------------------------------------------------

% caudal
peakAct_misMatch_caudal_prev_temp = squeeze(act_peak_prev_fin(:,:,1,:)).*(tType==1).*ifCaudal_used;
peakAct_misMatch_caudal_prev_temp_count = squeeze(sum(peakAct_misMatch_caudal_prev_temp,1));

peakAct_misMatch_caudal_prev_avg = peakAct_misMatch_caudal_prev_temp_count./tTypeCount_misMatch_caudal;

peakAct_match_caudal_prev_temp = squeeze(act_peak_prev_fin(:,:,1,:)).*(tType==2).*ifCaudal_used;
peakAct_match_caudal_prev_temp_count = squeeze(sum(peakAct_match_caudal_prev_temp,1));

peakAct_match_caudal_prev_avg = peakAct_match_caudal_prev_temp_count./tTypeCount_match_caudal;

% PRC
peakAct_misMatch_PRC_temp_prev = squeeze(act_peak_prev_fin(:,:,2,:)).*(tType==1).*ifPRC_used;
peakAct_misMatch_PRC_prev_temp_count = squeeze(sum(peakAct_misMatch_PRC_temp_prev,1));

peakAct_misMatch_PRC_prev_avg = peakAct_misMatch_PRC_prev_temp_count./tTypeCount_misMatch_PRC;

peakAct_match_PRC_prev_temp = squeeze(act_peak_prev_fin(:,:,2,:)).*(tType==2).*ifPRC_used;
peakAct_match_PRC__prev_temp_count = squeeze(sum(peakAct_match_PRC_prev_temp,1));

peakAct_match_PRC_prev_avg = peakAct_match_PRC__prev_temp_count./tTypeCount_match_PRC;

%--------------------------------------------------------------------------
%total
%--------------------------------------------------------------------------

% caudal
totalAct_misMatch_caudal_prev_temp = squeeze(act_total_prev_fin(:,:,1,:)).*(tType==1).*ifCaudal_used;
totalAct_misMatch_caudal_prev_temp_count = squeeze(sum(totalAct_misMatch_caudal_prev_temp,1));

totalAct_misMatch_caudal_prev_avg = totalAct_misMatch_caudal_prev_temp_count./tTypeCount_misMatch_caudal;

totalAct_match_caudal_prev_temp = squeeze(act_total_prev_fin(:,:,1,:)).*(tType==2).*ifCaudal_used;
totalAct_match_caudal_prev_temp_count = squeeze(sum(totalAct_match_caudal_prev_temp,1));

totalAct_match_caudal_prev_avg = totalAct_match_caudal_prev_temp_count./tTypeCount_match_caudal;

% PRC
totalAct_misMatch_PRC_prev_temp = squeeze(act_total_prev_fin(:,:,2,:)).*(tType==1).*ifPRC_used;
totalAct_misMatch_PRC_prev_temp_count = squeeze(sum(totalAct_misMatch_PRC_prev_temp,1));

totalAct_misMatch_PRC_prev_avg = totalAct_misMatch_PRC_prev_temp_count./tTypeCount_misMatch_PRC;

totalAct_match_PRC_prev_temp = squeeze(act_total_prev_fin(:,:,2,:)).*(tType==2).*ifPRC_used;
totalAct_match_PRC_prev_temp_count = squeeze(sum(totalAct_match_PRC_prev_temp,1));

totalAct_match_PRC_prev_avg = totalAct_match_PRC_prev_temp_count./tTypeCount_match_PRC;


%% familiarity difference that was used

familDiff_used_misMatch = familDiff_used.*(tType==1);
familDiff_used_match = familDiff_used.*(tType==2);

familDiff_used_misMatch_count = squeeze(sum(familDiff_used_misMatch,1));
familDiff_used_match_count = squeeze(sum(familDiff_used_match,1));

familDiff_used_misMatch_avg = familDiff_used_misMatch_count./tTypeCount_misMatch;
familDiff_used_match_avg = familDiff_used_match_count./tTypeCount_match;

%% whether PRC used

% previous stimuli
ifPRC_misMatch_temp_prev = squeeze(ifPRC_both(:,:,1,:)).*(tType==1);
ifPRC_misMatch_temp_prev_count = squeeze(sum(ifPRC_misMatch_temp_prev,1));

ifPRC_misMatch_prev_avg = ifPRC_misMatch_temp_prev_count./tTypeCount_misMatch;

ifPRC_match_temp_prev = squeeze(ifPRC_both(:,:,1,:)).*(tType==2);
ifPRC_match_temp_prev_count = squeeze(sum(ifPRC_match_temp_prev,1));

ifPRC_match_prev_avg = ifPRC_match_temp_prev_count./tTypeCount_match;

% new stimuli
ifPRC_misMatch_temp_new = squeeze(ifPRC_both(:,:,2,:)).*(tType==1);
ifPRC_misMatch_temp_new_count = squeeze(sum(ifPRC_misMatch_temp_new,1));

ifPRC_misMatch_new_avg = ifPRC_misMatch_temp_new_count./tTypeCount_misMatch;

ifPRC_match_temp_new = squeeze(ifPRC_both(:,:,2,:)).*(tType==2);
ifPRC_match_temp_new_count = squeeze(sum(ifPRC_match_temp_new,1));

ifPRC_match_new_avg = ifPRC_match_temp_new_count./tTypeCount_match;

% both stimuli
ifPRC_misMatch_temp_both = ifPRC_used.*(tType==1);
ifPRC_misMatch_temp_both_count = squeeze(sum(ifPRC_misMatch_temp_both,1));

ifPRC_misMatch_both_avg = ifPRC_misMatch_temp_both_count./tTypeCount_misMatch;

ifPRC_match_temp_both = ifPRC_used.*(tType==2);
ifPRC_match_temp_both_count = squeeze(sum(ifPRC_match_temp_both,1));

ifPRC_match_both_avg = ifPRC_match_temp_both_count./tTypeCount_match;



%% Selectivity Absolute

%**************************************************************************
% prev
%**************************************************************************

%--------------------------------------------------------------------------
% Caudal
%--------------------------------------------------------------------------

% misMatch
meanSelec_caudal_misMatch_prev = meanSelectivity_caudal_prev.*(tType==1).*ifCaudal_used;
meanSelec_caudal_misMatch_prev_count = squeeze(sum(meanSelec_caudal_misMatch_prev,1));

meanSelec_caudal_misMatch_prev_avg = meanSelec_caudal_misMatch_prev_count./tTypeCount_misMatch_caudal;

% match
meanSelec_caudal_match_prev = meanSelectivity_caudal_prev.*(tType==2).*ifCaudal_used;
meanSelec_caudal_match_prev_count = squeeze(sum(meanSelec_caudal_match_prev,1));

meanSelec_caudal_match_prev_avg = meanSelec_caudal_match_prev_count./tTypeCount_match_caudal;
%--------------------------------------------------------------------------
% PRC
%--------------------------------------------------------------------------

% misMatch
meanSelec_PRC_misMatch_prev = meanSelectivity_PRC_prev.*(tType==1).*ifPRC_used;
meanSelec_PRC_misMatch_prev_count = squeeze(sum(meanSelec_PRC_misMatch_prev,1));

meanSelec_PRC_misMatch_prev_avg = meanSelec_PRC_misMatch_prev_count./tTypeCount_misMatch_PRC;

% match
meanSelec_PRC_match_prev = meanSelectivity_PRC_prev.*(tType==2).*ifPRC_used;
meanSelec_PRC_match_prev_count = squeeze(sum(meanSelec_PRC_match_prev,1));

meanSelec_PRC_match_prev_avg = meanSelec_PRC_match_prev_count./tTypeCount_match_PRC;

%**************************************************************************
% new
%**************************************************************************

%--------------------------------------------------------------------------
% Caudal
%--------------------------------------------------------------------------

% misMatch
meanSelec_caudal_misMatch_new = meanSelectivity_caudal_new.*(tType==1).*ifCaudal_used;
meanSelec_caudal_misMatch_new_count = squeeze(sum(meanSelec_caudal_misMatch_new,1));

meanSelec_caudal_misMatch_new_avg = meanSelec_caudal_misMatch_new_count./tTypeCount_misMatch_caudal;

% match
meanSelec_caudal_match_new = meanSelectivity_caudal_new.*(tType==2).*ifCaudal_used;
meanSelec_caudal_match_new_count = squeeze(sum(meanSelec_caudal_match_new,1));

meanSelec_caudal_match_new_avg = meanSelec_caudal_match_new_count./tTypeCount_match_caudal;
%--------------------------------------------------------------------------
% PRC
%--------------------------------------------------------------------------

% misMatch
meanSelec_PRC_misMatch_new = meanSelectivity_PRC_new.*(tType==1).*ifPRC_used;
meanSelec_PRC_misMatch_new_count = squeeze(sum(meanSelec_PRC_misMatch_new,1));

meanSelec_PRC_misMatch_new_avg = meanSelec_PRC_misMatch_new_count./tTypeCount_misMatch_PRC;

% match
meanSelec_PRC_match_new = meanSelectivity_PRC_new.*(tType==2).*ifPRC_used;
meanSelec_PRC_match_new_count = squeeze(sum(meanSelec_PRC_match_new,1));

meanSelec_PRC_match_new_avg = meanSelec_PRC_match_new_count./tTypeCount_match_PRC;

%% Famil Differences, actual

familDifferences_misMatch_caudal_temp = squeeze(sum(familDifferences_misMatch_caudal,1));
familDiffs_misMatch_caudal_mean = familDifferences_misMatch_caudal_temp./trial_misMatch_caudal;

familDifferences_match_caudal_temp = squeeze(sum(familDifferences_match_caudal,1));
familDiffs_match_caudal_mean = familDifferences_match_caudal_temp./trial_match_caudal;

familDifferences_misMatch_PRC_temp = squeeze(sum(familDifferences_misMatch_PRC,1));
familDiffs_misMatch_PRC_mean = familDifferences_misMatch_PRC_temp(:,:)./trial_misMatch_PRC(:,:);

familDifferences_match_PRC_temp = squeeze(sum(familDifferences_match_PRC,1));
familDiffs_match_PRC_mean = familDifferences_match_PRC_temp(:,:)./trial_match_PRC(:,:);

%% calculating total fixations

% lesion and control by session
fixations_misMatch_temp = squeeze(sum(fixations_total.*(tType==1),1));
fixations_misMatch_mean = fixations_misMatch_temp(:,:)./tTypeCount_misMatch(:,:);

fixations_match_temp = squeeze(sum(fixations_total.*(tType==2),1));
fixations_match_mean = fixations_match_temp(:,:)./tTypeCount_match(:,:);

% divide control by PRC and caudal

fixations_caudal_misMatch_temp = squeeze(sum(fixations_total.*(tType==1).*ifCaudal_used,1));
fixations_caudal_misMatch_mean = fixations_caudal_misMatch_temp(:,:)./trial_misMatch_caudal(:,:);

fixations_caudal_match_temp = squeeze(sum(fixations_total.*(tType==2).*ifCaudal_used,1));
fixations_caudal_match_mean = fixations_caudal_match_temp(:,:)./trial_match_caudal(:,:);

fixations_PRC_misMatch_temp = squeeze(sum(fixations_total.*(tType==1).*ifPRC_used,1));
fixations_PRC_misMatch_mean = fixations_PRC_misMatch_temp(:,:)./trial_misMatch_PRC(:,:);

fixations_PRC_match_temp = squeeze(sum(fixations_total.*(tType==2).*ifPRC_used,1));
fixations_PRC_match_mean = fixations_PRC_match_temp(:,:)./trial_match_PRC(:,:);


%% STD of familDifferences

familDiffs_misMatch_caudal_var_toSqr = zeros(nSess,p.nTrials);
familDiffs_match_caudal_var_toSqr = zeros(nSess,p.nTrials);
familDiffs_misMatch_PRC_var_toSqr = zeros(nSess,p.nTrials);
familDiffs_match_PRC_var_toSqr = zeros(nSess,p.nTrials);



for session = 1:nSess
    for trial = 1:p.nTrials
        familDiffs_misMatch_caudal_std_temp = 0;
        familDiffs_match_caudal_std_temp = 0;
        familDiffs_misMatch_PRC_std_temp = 0;
        familDiffs_match_PRC_std_temp = 0;
        for rat = 1:numRats
            
            % find sum of squared differences from the mean
            if trial_misMatch_caudal_if(rat,session,trial)
                familDiffs_misMatch_caudal_std_temp = familDiffs_misMatch_caudal_std_temp + (familDifferences_misMatch_caudal(rat,session,trial) - familDiffs_misMatch_caudal_mean(session,trial)).^2;
            end
            
            if trial_match_caudal_if(rat,session,trial)
                familDiffs_match_caudal_std_temp = familDiffs_match_caudal_std_temp + (familDifferences_match_caudal(rat,session,trial) - familDiffs_match_caudal_mean(session,trial)).^2;
            end
            
            if trial_misMatch_PRC_if(rat,session,trial)
                familDiffs_misMatch_PRC_std_temp = familDiffs_misMatch_PRC_std_temp + (familDifferences_misMatch_PRC(rat,session,trial) - familDiffs_misMatch_PRC_mean(session,trial)).^2;
            end
            
            if trial_match_PRC_if(rat,session,trial)
                familDiffs_match_PRC_std_temp = familDiffs_match_PRC_std_temp + (familDifferences_match_PRC(rat,session,trial) - familDiffs_match_PRC_mean(session,trial)).^2;
            end
            
        end
        
        % divide by n-1
        familDiffs_misMatch_caudal_var_toSqr(session,trial) = familDiffs_misMatch_caudal_std_temp / (trial_misMatch_caudal(session,trial) -1);
        familDiffs_match_caudal_var_toSqr(session,trial) = familDiffs_match_caudal_std_temp / (trial_match_caudal(session,trial) -1);
        familDiffs_misMatch_PRC_var_toSqr(session,trial) = familDiffs_misMatch_PRC_std_temp / (trial_misMatch_PRC(session,trial) -1);
        familDiffs_match_PRC_var_toSqr(session,trial) = familDiffs_match_PRC_std_temp / (trial_match_PRC(session,trial) -1);
    end
end


% put in csv
% fName = 'results.csv';
% fileToWrite = fopen(fName,'w');
%
% fprintf()
%


% squareroot those variances to get the standard deviation
familDiffs_misMatch_caudal_std = sqrt(familDiffs_misMatch_caudal_var_toSqr);
familDiffs_match_caudal_std = sqrt(familDiffs_match_caudal_var_toSqr);
familDiffs_misMatch_PRC_std = sqrt(familDiffs_misMatch_PRC_var_toSqr);
familDiffs_match_PRC_std = sqrt(familDiffs_match_PRC_var_toSqr);

familDiffs_misMatch_caudal_sem = familDiffs_misMatch_caudal_std ./ sqrt(trial_misMatch_caudal(:,:));
familDiffs_match_caudal_sem = familDiffs_match_caudal_std ./ sqrt(trial_match_caudal(:,:));
familDiffs_misMatch_PRC_sem = familDiffs_misMatch_PRC_std ./ sqrt(trial_misMatch_PRC(:,:));
familDiffs_match_PRC_sem = familDiffs_match_PRC_std ./ sqrt(trial_match_PRC(:,:));

trials = 1:p.nTrials;
familDiffs_match_caudal_err = [familDiffs_match_caudal_mean + familDiffs_match_caudal_sem; familDiffs_match_caudal_mean - familDiffs_match_caudal_sem];
familDiffs_misMatch_caudal_err = [familDiffs_misMatch_caudal_mean + familDiffs_misMatch_caudal_sem; familDiffs_misMatch_caudal_mean - familDiffs_misMatch_caudal_sem];
familDiffs_match_PRC_err = [familDiffs_match_PRC_mean + familDiffs_match_PRC_sem; familDiffs_match_PRC_mean - familDiffs_match_PRC_sem];
familDiffs_misMatch_PRC_err = [familDiffs_misMatch_PRC_mean + familDiffs_misMatch_PRC_sem; familDiffs_misMatch_PRC_mean - familDiffs_misMatch_PRC_sem];


%% looking at more raw hit and FA rats


dPrime_first_raw = norminv(hitRate_first) - norminv(FARate_first);
dPrime_second_raw = norminv(hitRate_second) - norminv(FARate_second);

dPrime_first_raw_mean = mean(dPrime_first_raw);
dPrime_second_raw_mean= mean(dPrime_second_raw);

first_std = std(dPrime_first_raw);
first_sem = first_std/sqrt(numRats);

second_std = std(dPrime_second_raw);
second_sem = second_std/sqrt(numRats);

dPrime_first_err = first_sem;
dPrime_second_err = second_sem;

if any(isnan(dPrime_first_err))
    dPrime_first_err(:) = 0;
end

if any(isnan(dPrime_second_err))
    dPrime_second_err(:) = 0;
end


dPrime_first_adj_all = norminv(hitRate_first_adj_all) - norminv(FARate_first_adj_all);
dPrime_second_adj_all = norminv(hitRate_second_adj_all) - norminv(FARate_second_adj_all);

dPrime_first_adj_all_mean = mean(dPrime_first_adj_all);
dPrime_second_adj_all_mean = mean(dPrime_second_adj_all);

first_adj_all_std = std(dPrime_first_adj_all);
first_adj_all_sem = first_adj_all_std/sqrt(numRats);

second_adj_all_std = std(dPrime_second_adj_all);
second_adj_all_sem = second_adj_all_std/sqrt(numRats);

dPrime_first_adj_all_err = first_adj_all_sem;
dPrime_second_adj_all_err = second_adj_all_sem;


dPrime_first_adj_some = norminv(hitRate_first_adj_some) - norminv(FARate_first_adj_some);
dPrime_second_adj_some = norminv(hitRate_second_adj_some) - norminv(FARate_second_adj_some);

dPrime_first_adj_some_mean = mean(dPrime_first_adj_some);
dPrime_second_adj_some_mean = mean(dPrime_second_adj_some);

first_adj_some_std = std(dPrime_first_adj_some);
first_adj_some_sem = first_adj_some_std/sqrt(numRats);

second_adj_some_std = std(dPrime_second_adj_some);
second_adj_some_sem = second_adj_some_std/sqrt(numRats);

dPrime_first_adj_some_err = first_adj_some_sem;
dPrime_second_adj_some_err = second_adj_some_sem;

%% dPrime of total sessions

syms dp
dP_adj_all = zeros(size(hitRate_all));
for i = 1:size(hitRate_all,1)
    for j = 1:size(hitRate_all,2)
        H = hitRate_all(i,j);
        k = sqrt(2).*norminv(FARate_all(i,j)/2);
        eqn = H == normcdf((k+dp)/sqrt(2)) + normcdf((k-dp)/sqrt(2));
        
        dP_adj_all(i,j) = abs(solve(eqn,dp));
    end
end

% dP_adj_all = norminv(hitRate_all) - norminv(FARate_all);

dP_adj_all_mean = mean(dP_adj_all,1);

dP_adj_all_std = std(dP_adj_all);
dP_adj_all_sem = dP_adj_all_std/sqrt(numRats);

% d' of just the comparison trials
syms dp
dP_all_comp = zeros(size(hitRate_all_comp));
for i = 1:size(hitRate_all_comp,1)
    for j = 1:size(hitRate_all_comp,2)
        H = hitRate_all_comp(i,j);
        k = sqrt(2).*norminv(FARate_all_comp(i,j)/2);
        eqn = H == normcdf((k+dp)/sqrt(2)) + normcdf((k-dp)/sqrt(2));
        
        dP_all_comp(i,j) = abs(solve(eqn,dp));
    end
end

% dP_all_comp = norminv(hitRate_all_comp) - norminv(FARate_all_comp);

dP_all_comp_mean = mean(dP_all_comp,1);

dP_all_comp_std = std(dP_all_comp);
dP_all_comp_sem = dP_all_comp_std/sqrt(numRats);

%% Most important graphs!
close all;
%% famil difference
if ~onlyFigure
    
    
    %--------------------------------------------------------------------------
    % define axis properties
    %--------------------------------------------------------------------------
    maxY = max(max([familDiffs_misMatch_caudal_mean, familDiffs_match_caudal_mean, familDiffs_match_PRC_mean, familDiffs_misMatch_PRC_mean]));
    minY_familDiffs = 0;
    maxY_familDiffs = maxY;
    axisBounds_familDiffs = [0,p.nTrials,minY_familDiffs,maxY_familDiffs];
    
    %--------------------------------------------------------------------------
    % First, make plots for lesion sessions
    %--------------------------------------------------------------------------
    
    % CAUDAL graphs
    figs(1) = figure;
    hold on
    plot(familDiffs_match_caudal_mean(1,:), 'color', 'g')
    plot(familDiffs_match_caudal_mean(2,:), 'color', 'r')
    plot(familDiffs_match_caudal_mean(3,:), 'color', 'b')
    plot(familDiffs_misMatch_caudal_mean(1,:), 'g', 'LineStyle', '--')
    plot(familDiffs_misMatch_caudal_mean(2,:), 'color','r', 'LineStyle', '--')
    plot(familDiffs_misMatch_caudal_mean(3,:), 'color','b', 'LineStyle', '--')
    
    plotshaded(trials,[familDiffs_misMatch_caudal_err(1,:);familDiffs_misMatch_caudal_err(7,:)],'g')
    plotshaded(trials,[familDiffs_misMatch_caudal_err(2,:);familDiffs_misMatch_caudal_err(8,:)],'r')
    plotshaded(trials,[familDiffs_misMatch_caudal_err(3,:);familDiffs_misMatch_caudal_err(9,:)],'b')
    
    plotshaded(trials,[familDiffs_match_caudal_err(1,:);familDiffs_match_caudal_err(7,:)],'g')
    plotshaded(trials,[familDiffs_match_caudal_err(2,:);familDiffs_match_caudal_err(8,:)],'r')
    plotshaded(trials,[familDiffs_match_caudal_err(3,:);familDiffs_match_caudal_err(9,:)],'b')
    
    xlabel('trial');
    ylabel('familiarity difference');
    legend('low1,match', 'high,match', 'low2,match', 'low1,misMatch', 'high,misMatch', 'low2,misMatch','Location','best');
    title({'Familiarity diffs, lesion'});
    axis(axisBounds_familDiffs);
    
    % NO PRC
    
    saveas(figs(1),[saveFolder, '/famDiff_les'], 'fig');
    saveas(figs(1),[saveFolder, '/famDiff_les'], 'jpg');
    
    
    %--------------------------------------------------------------------------
    % next, plot CONTROL sessions
    %--------------------------------------------------------------------------
    
    % caudal
    figs(2) = figure;
    subplot(1,2,1)
    hold on
    plot(familDiffs_match_caudal_mean(4,:), 'color', 'g')
    plot(familDiffs_match_caudal_mean(5,:), 'color', 'r')
    plot(familDiffs_match_caudal_mean(6,:), 'color', 'b')
    plot(familDiffs_misMatch_caudal_mean(4,:), 'g', 'LineStyle', '--')
    plot(familDiffs_misMatch_caudal_mean(5,:), 'color','r', 'LineStyle', '--')
    plot(familDiffs_misMatch_caudal_mean(6,:), 'color','b', 'LineStyle', '--')
    
    
    plotshaded(trials,[familDiffs_misMatch_caudal_err(4,:);familDiffs_misMatch_caudal_err(10,:)],'g')
    plotshaded(trials,[familDiffs_misMatch_caudal_err(5,:);familDiffs_misMatch_caudal_err(11,:)],'r')
    plotshaded(trials,[familDiffs_misMatch_caudal_err(6,:);familDiffs_misMatch_caudal_err(12,:)],'b')
    
    plotshaded(trials,[familDiffs_match_caudal_err(4,:);familDiffs_match_caudal_err(10,:)],'g')
    plotshaded(trials,[familDiffs_match_caudal_err(5,:);familDiffs_match_caudal_err(11,:)],'r')
    plotshaded(trials,[familDiffs_match_caudal_err(6,:);familDiffs_match_caudal_err(12,:)],'b')
    
    xlabel('trial');
    ylabel('familiarity difference');
    legend('low1,match', 'high,match', 'low2,match', 'low1,misMatch', 'high,misMatch', 'low2,misMatch','Location','best');
    title({'Familiarity diffs, CONTROL, CAUDAL'});
    axis(axisBounds_familDiffs);
    
    %--------------------------------------------------------------------------
    % PRC
    subplot(1,2,2)
    hold on
    plot(familDiffs_match_PRC_mean(4,:), 'color', 'g')
    plot(familDiffs_match_PRC_mean(5,:), 'color', 'r')
    plot(familDiffs_match_PRC_mean(6,:), 'color', 'b')
    
    plot(familDiffs_misMatch_PRC_mean(4,:), 'g', 'LineStyle', '--')
    plot(familDiffs_misMatch_PRC_mean(5,:), 'color','r', 'LineStyle', '--')
    plot(familDiffs_misMatch_PRC_mean(6,:), 'color','b', 'LineStyle', '--')
    
    plotshaded(trials,[familDiffs_misMatch_PRC_err(4,:);familDiffs_misMatch_PRC_err(10,:)],'g')
    plotshaded(trials,[familDiffs_misMatch_PRC_err(5,:);familDiffs_misMatch_PRC_err(11,:)],'r')
    plotshaded(trials,[familDiffs_misMatch_PRC_err(6,:);familDiffs_misMatch_PRC_err(12,:)],'r')
    
    plotshaded(trials,[familDiffs_match_PRC_err(4,:);familDiffs_match_PRC_err(10,:)],'g')
    plotshaded(trials,[familDiffs_match_PRC_err(5,:);familDiffs_match_PRC_err(11,:)],'r')
    plotshaded(trials,[familDiffs_match_PRC_err(6,:);familDiffs_match_PRC_err(12,:)],'r')
    
    xlabel('trial');
    ylabel('familiarity difference');
    legend('low1,match', 'high,match', 'low2,match', 'low1,misMatch', 'high,misMatch', 'low2,misMatch','Location','best');
    title({'Familiarity diffs, CONTROL, PRC'});
    axis(axisBounds_familDiffs);
    
    saveas(figs(2),[saveFolder, '/famDiff_cntr'], 'fig');
    saveas(figs(2),[saveFolder, '/famDiff_cntr'], 'jpg');
    
    
    %% d' from raw scores
    
    %--------------------------------------------------------------------------
    % define axis properties
    %--------------------------------------------------------------------------
    % minY_dPrime = 0;
    % maxY_dPrime = 1.8;
    % axisBounds_familDiffs = [0,p.nTrials,minY_dPrime,maxY_dPrime];
    %
    % %
    % figs(3) = figure;
    % subplot(1,2,1)
    % barweb([dPrime_first_raw_mean(2), dPrime_second_raw_mean(2) ; dPrime_first_raw_mean(1), dPrime_second_raw_mean(1)], [dPrime_first_err(2), dPrime_second_err(2) ; dPrime_first_err(1), dPrime_second_err(1)], [], {'high', 'low'})
    % xlabel('stim ambiguity');
    % ylabel('d''');
    % legend('first Half', 'Second Half')
    % title({'dPrime of Lesion'})
    %
    % subplot(1,2,2)
    % barweb([dPrime_first_raw_mean(4), dPrime_second_raw_mean(4) ; dPrime_first_raw_mean(3), dPrime_second_raw_mean(3)], [dPrime_first_err(4), dPrime_second_err(4) ; dPrime_first_err(3), dPrime_second_err(3)], [], {'high', 'low'})
    % xlabel('stim ambiguity');
    % ylabel('d''');
    % legend('first Half', 'Second Half')
    % title({'dPrime of Control'})
    %
    % saveas(figs(3),[saveFolder, '/dPrime'], 'fig');
    % saveas(figs(3),[saveFolder, '/dPrime'], 'jpg');
    %
    
    % %% aboslute famil
    %
    % %--------------------------------------------------------------------------
    % % define axis properties
    % %--------------------------------------------------------------------------
    % minY_abslFamil = 0;
    % maxY_abslFamil = .004;
    % axisBounds_abslFamil = [0,p.nTrials,minY_abslFamil,maxY_abslFamil];
    %
    %
    % %--------------------------------------------------------------------------
    % % LESION
    % %--------------------------------------------------------------------------
    %
    % % caudal, prev
    % figs(4) = figure;
    % subplot(1,2,1)
    % hold on
    % plot(meanSelec_caudal_match_prev_avg(1,:), 'color', 'g')
    % plot(meanSelec_caudal_match_prev_avg(2,:), 'color', 'r')
    % plot(meanSelec_caudal_misMatch_prev_avg(1,:), 'g', 'LineStyle', '--')
    % plot(meanSelec_caudal_misMatch_prev_avg(2,:), 'color','r', 'LineStyle', '--')
    % legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
    % title({'Absolute Familiarity, caudal, PREV, LESION'});
    % axis(axisBounds_abslFamil)
    %
    % % caudal, new
    % subplot(1,2,2); hold on
    % plot(meanSelec_caudal_match_new_avg(1,:), 'color', 'g')
    % plot(meanSelec_caudal_match_prev_avg(2,:), 'color', 'r')
    % plot(meanSelec_caudal_misMatch_prev_avg(1,:), 'g', 'LineStyle', '--')
    % plot(meanSelec_caudal_misMatch_prev_avg(2,:), 'color','r', 'LineStyle', '--')
    % legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
    % title({'Absolute Familiarity,NEW, caudal, LESION'});
    % axis(axisBounds_abslFamil)
    %
    % saveas(figs(4),[saveFolder, '/abFam_caud_les'], 'fig');
    % saveas(figs(4),[saveFolder, '/abFam_caud_les'], 'jpg');
    %
    % % no PRC to plot
    %
    % %--------------------------------------------------------------------------
    % % next, plot CONTROL sessions
    % %--------------------------------------------------------------------------
    %
    % % caudal, prev
    % figs(5) = figure;
    % subplot(1,2,1)
    % hold on
    % plot(meanSelec_caudal_match_prev_avg(3,:), 'color', 'g')
    % plot(meanSelec_caudal_match_prev_avg(4,:), 'color', 'r')
    % plot(meanSelec_caudal_misMatch_prev_avg(3,:), 'g', 'LineStyle', '--')
    % plot(meanSelec_caudal_misMatch_prev_avg(4,:), 'color','r', 'LineStyle', '--')
    % legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
    % title({'Absolute Familiarity, PREV, caudal, CONTROl'});
    % axis(axisBounds_abslFamil)
    %
    % % caudal, new
    % subplot(1,2,2); hold on
    % plot(meanSelec_caudal_match_new_avg(3,:), 'color', 'g')
    % plot(meanSelec_caudal_match_new_avg(4,:), 'color', 'r')
    % plot(meanSelec_caudal_misMatch_new_avg(3,:), 'g', 'LineStyle', '--')
    % plot(meanSelec_caudal_misMatch_new_avg(4,:), 'color','r', 'LineStyle', '--')
    % legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
    % title({'Absolute Familiarity,NEW, caudal, CONTROL'});
    % axis(axisBounds_abslFamil)
    %
    % saveas(figs(5),[saveFolder, '/bFam_cud_cntr'], 'fig');
    % saveas(figs(5),[saveFolder, '/bFam_cud_cntr'], 'jpg');
    %
    % % PRC
    % figs(6) = figure;
    % subplot(1,2,1)
    % hold on
    % plot(meanSelec_PRC_match_prev_avg(3,:), 'color', 'g')
    % plot(meanSelec_PRC_match_prev_avg(4,:), 'color', 'r')
    % plot(meanSelec_PRC_misMatch_prev_avg(3,:), 'g', 'LineStyle', '--')
    % plot(meanSelec_PRC_misMatch_prev_avg(4,:), 'color','r', 'LineStyle', '--')
    % legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
    % title({'Absolute Familiarity, PREV, PRC'});
    % axis(axisBounds_abslFamil)
    %
    % subplot(1,2,2); hold on
    % plot(meanSelec_PRC_match_new_avg(3,:), 'color', 'g')
    % plot(meanSelec_PRC_match_new_avg(4,:), 'color', 'r')
    % plot(meanSelec_PRC_misMatch_new_avg(3,:), 'g', 'LineStyle', '--')
    % plot(meanSelec_PRC_misMatch_new_avg(4,:), 'color','r', 'LineStyle', '--')
    % legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
    % title({'Absolute Familiarity,NEW, PRC'});
    % axis(axisBounds_abslFamil)
    %
    % saveas(figs(6),[saveFolder, '/abFam_PRC'], 'fig');
    % saveas(figs(6),[saveFolder, '/abFam_PRC'], 'jpg');
    %
    %
    % %% feats sampled per comparison
    %
    % figs(7) = figure; hold on
    % plot(fixByComp_prev(1,:,:), 'color', 'g')
    % plot(fixByComp_prev(2,:,:), 'color', 'r')
    % plot(fixByComp_prev(3,:,:), 'g', 'LineStyle', '--')
    % plot(fixByComp_prev(4,:,:), 'color','r', 'LineStyle', '--')
    % xlabel('trial');
    % ylabel('features sampled');
    % legend('lesion,low', 'lesion,high', 'control,low', 'control,high');
    % title({'Avg Features per Comparison, PREV'});
    %
    % % plot(squeeze(fixByComp_prev_mean_plot(1,2,:)))
    % saveas(figs(7),[saveFolder, '/nFixPerComp'], 'fig');
    % saveas(figs(7),[saveFolder, '/nFixPerComp'], 'jpg');
    %
    %
    % figs(8) = figure; hold on
    % plot(featSamplePerTrial_mean(1,:), 'color', 'g')
    % plot(featSamplePerTrial_mean(2,:), 'color', 'g', 'LineStyle', '--')
    % plot(featSamplePerTrial_mean(3,:), 'r')
    % plot(featSamplePerTrial_mean(4,:), 'color','r', 'LineStyle', '--')
    % xlabel('trial');
    % ylabel('total features sampled');
    % legend('lesion,low', 'lesion,high', 'control,low', 'control,high');
    % title({'total features sampled at each trial'});
    %
    % saveas(figs(8),[saveFolder, '/nTFtSmpld'], 'fig');
    % saveas(figs(8),[saveFolder, '/nTFtSmpld'], 'jpg');
    %
    %
    % %%
    % figs(9) = figure;
    % subplot(1,2,1)
    % hold on
    % plot(featSamplePerTrial_match_mean(1,:), 'color', 'g')
    % plot(featSamplePerTrial_match_mean(2,:), 'color', 'r')
    % plot(featSamplePerTrial_misMatch_mean(1,:), 'g', 'LineStyle', '--')
    % plot(featSamplePerTrial_misMatch_mean(2,:), 'color','r', 'LineStyle', '--')
    % xlabel('trial');
    % ylabel('total features sampled');
    % legend('match,low', 'match,high', 'misMatch,low', 'misMatch,high');
    % title({'total features sampled at each trial, LESION'});
    %
    % subplot(1,2,2)
    % hold on
    % plot(featSamplePerTrial_match_mean(3,:), 'color', 'g')
    % plot(featSamplePerTrial_match_mean(4,:), 'color', 'r')
    % plot(featSamplePerTrial_misMatch_mean(3,:), 'g', 'LineStyle', '--')
    % plot(featSamplePerTrial_misMatch_mean(4,:), 'color','r', 'LineStyle', '--')
    % xlabel('trial');
    % ylabel('total features sampled');
    % legend('match,low', 'match,high', 'misMatch,low', 'misMatch,high');
    % title({'total features sampled at each trial, CONTROL'});
    %
    % saveas(figs(9),[saveFolder, '/nFtsSmpld_Trl'], 'fig');
    % saveas(figs(9),[saveFolder, '/nFtsSmpld_Trl'], 'jpg');
    %
    %
    % %% total fixations
    % minY_fixations = 0;
    % maxY_fixations = 60;
    % axisBounds_fixations = [0,p.nTrials,minY_fixations,maxY_fixations];
    %
    %
    % figs(10) = figure;
    % subplot(1,2,1)
    % hold on
    % plot(fixations_match_mean(1,:), 'color', 'g')
    % plot(fixations_match_mean(2,:), 'color', 'r')
    % plot(fixations_misMatch_mean(1,:), 'g', 'LineStyle', '--')
    % plot(fixations_misMatch_mean(2,:), 'color','r', 'LineStyle', '--')
    % xlabel('trial');
    % ylabel('total features sampled');
    % legend('match,low', 'match,high', 'misMatch,low', 'misMatch,high');
    % title({'total fixtions by trial, LESION'});
    % axis(axisBounds_fixations)
    %
    %
    % subplot(1,2,2)
    % hold on
    % plot(fixations_match_mean(3,:), 'color', 'g')
    % plot(fixations_match_mean(4,:), 'color', 'r')
    % plot(fixations_misMatch_mean(3,:), 'g', 'LineStyle', '--')
    % plot(fixations_misMatch_mean(4,:), 'color','r', 'LineStyle', '--')
    % xlabel('trial');
    % ylabel('total features sampled');
    % legend('match,low', 'match,high', 'misMatch,low', 'misMatch,high');
    % title({'total fixtions by trial, CONTROL'});
    % axis(axisBounds_fixations)
    %
    % saveas(figs(10),[saveFolder, '/nFixns'], 'fig');
    % saveas(figs(10),[saveFolder, '/nFixns'], 'jpg');
    %
    %
    % % divide control by PRC and Caudal
    % figs(11) = figure;
    % subplot(1,2,1)
    % hold on
    % plot(fixations_caudal_match_mean(3,:), 'color', 'g')
    % plot(fixations_caudal_match_mean(4,:), 'color', 'r')
    % plot(fixations_caudal_misMatch_mean(3,:), 'g', 'LineStyle', '--')
    % plot(fixations_caudal_misMatch_mean(4,:), 'color','r', 'LineStyle', '--')
    % xlabel('trial');
    % ylabel('total features sampled');
    % legend('match,low', 'match,high', 'misMatch,low', 'misMatch,high');
    % title({'total fixtions by trial, caudal'});
    % axis(axisBounds_fixations)
    %
    %
    % subplot(1,2,2)
    % hold on
    % plot(fixations_PRC_match_mean(3,:), 'color', 'g')
    % plot(fixations_PRC_match_mean(4,:), 'color', 'r')
    % plot(fixations_PRC_misMatch_mean(3,:), 'g', 'LineStyle', '--')
    % plot(fixations_PRC_misMatch_mean(4,:), 'color','r', 'LineStyle', '--')
    % xlabel('trial');
    % ylabel('total features sampled');
    % legend('match,low', 'match,high', 'misMatch,low', 'misMatch,high');
    % title({'total fixtions by trial, PRC'});
    % axis(axisBounds_fixations)
    %
    % saveas(figs(11),[saveFolder, '/nFixns_cntrl'], 'fig');
    % saveas(figs(11),[saveFolder, '/nFixns_cntrl'], 'jpg');
    %
    % %% absolFamil, better for looking at current vs new
    %
    % %--------------------------------------------------------------------------
    % % define axis properties
    % %--------------------------------------------------------------------------
    % minY_abslFamil_newPrev = minY_abslFamil;
    % maxY_abslFamil_newPrev = maxY_abslFamil;
    % axisBounds_abslFamil_newPrev = [0,p.nTrials,minY_abslFamil_newPrev,maxY_abslFamil_newPrev];
    %
    % %--------------------------------------------------------------------------
    % % LESION
    % %--------------------------------------------------------------------------
    %
    % % caudal
    % figs(12) = figure;
    % subplot(1,2,1)
    % hold on
    % plot(meanSelec_caudal_match_prev_avg(1,:), 'color', 'g')
    % plot(meanSelec_caudal_match_new_avg(1,:), 'color', 'g', 'LineStyle', '--');
    % plot(meanSelec_caudal_misMatch_prev_avg(1,:), 'r');
    % plot(meanSelec_caudal_misMatch_new_avg(1,:), 'r', 'LineStyle', '--');
    % legend('match,prev', 'match,new', 'misMatch,prev', 'misMatch,new');
    % title({'Absolute Famil, LA, LESION, caudal'});
    % axis(axisBounds_abslFamil_newPrev);
    %
    % subplot(1,2,2)
    % hold on
    % plot(meanSelec_caudal_match_prev_avg(2,:), 'color', 'g')
    % plot(meanSelec_caudal_match_new_avg(2,:), 'color', 'g', 'LineStyle', '--');
    % plot(meanSelec_caudal_misMatch_prev_avg(2,:), 'r');
    % plot(meanSelec_caudal_misMatch_new_avg(2,:), 'r', 'LineStyle', '--');
    % legend('match,prev', 'match,new', 'misMatch,prev', 'misMatch,new');
    % title({'Absolute Familiarity, HA, LESION, caudal'});
    % axis(axisBounds_abslFamil_newPrev);
    %
    % saveas(figs(12),[saveFolder, '/abFam_caud_NEW'], 'fig');
    % saveas(figs(12),[saveFolder, '/abFam_caud_NEW'], 'jpg');
    %
    % % no PRC
    %
    % %--------------------------------------------------------------------------
    % % CONTROL
    % %--------------------------------------------------------------------------
    %
    % % caudal
    % figs(13) = figure;
    % subplot(1,2,1)
    % hold on
    % plot(meanSelec_caudal_match_prev_avg(3,:), 'color', 'g')
    % plot(meanSelec_caudal_match_new_avg(3,:), 'color', 'g', 'LineStyle', '--');
    % plot(meanSelec_caudal_misMatch_prev_avg(3,:), 'r');
    % plot(meanSelec_caudal_misMatch_new_avg(3,:), 'r', 'LineStyle', '--');
    % legend('match,prev', 'match,new', 'misMatch,prev', 'misMatch,new');
    % title({'Absolute Famil, LA, CONTROL, caudal'});
    % axis(axisBounds_abslFamil_newPrev);
    %
    % subplot(1,2,2)
    % hold on
    % plot(meanSelec_caudal_match_prev_avg(4,:), 'color', 'g')
    % plot(meanSelec_caudal_match_new_avg(4,:), 'color', 'g', 'LineStyle', '--');
    % plot(meanSelec_caudal_misMatch_prev_avg(4,:), 'r');
    % plot(meanSelec_caudal_misMatch_new_avg(4,:), 'r', 'LineStyle', '--');
    % legend('match,prev', 'match,new', 'misMatch,prev', 'misMatch,new');
    % title({'Absolute Familiarity, HA, CONTROL, caudal'});
    % axis(axisBounds_abslFamil_newPrev);
    %
    % saveas(figs(13),[saveFolder, '/abFam_cud_C_N'], 'fig');
    % saveas(figs(13),[saveFolder, '/abFam_cud_C_N'], 'jpg');
    %
    %
    % % PRC
    % figs(14) = figure;
    % subplot(1,2,1)
    % hold on
    % plot(meanSelec_PRC_match_prev_avg(3,:), 'color', 'g')
    % plot(meanSelec_PRC_match_new_avg(3,:), 'color', 'g', 'LineStyle', '--');
    % plot(meanSelec_PRC_misMatch_prev_avg(3,:), 'r');
    % plot(meanSelec_PRC_misMatch_new_avg(3,:), 'r', 'LineStyle', '--');
    % legend('match,prev', 'match,new', 'misMatch,prev', 'misMatch,new');
    % title({'Absolute Familiarity, LA, CONTROL, PRC'});
    % axis(axisBounds_abslFamil_newPrev);
    %
    % subplot(1,2,2)
    % hold on
    % plot(meanSelec_PRC_match_prev_avg(4,:), 'color', 'g')
    % plot(meanSelec_PRC_match_new_avg(4,:), 'color', 'g', 'LineStyle', '--');
    % plot(meanSelec_PRC_misMatch_prev_avg(4,:), 'r');
    % plot(meanSelec_PRC_misMatch_new_avg(4,:), 'r', 'LineStyle', '--');
    % legend('match,prev', 'match,new', 'misMatch,prev', 'misMatch,new');
    % title({'Absolute Familiarity, HA, CONTROL, PRC'});
    % axis(axisBounds_abslFamil_newPrev);
    %
    % saveas(figs(14),[saveFolder, '/abFam_PRC_NEW'], 'fig');
    % saveas(figs(14),[saveFolder, '/abFam_PRC_NEW'], 'jpg');
    
    %% peak by match/misMatch
    
    %--------------------------------------------------------------------------
    % define axis properties
    %--------------------------------------------------------------------------
    minY_peak_tType_caudal = 0;
    maxY_peak_tType_caudal = p.sizeOfPeak;
    axisBounds_peak_tType_caudal = [0,p.nTrials,minY_peak_tType_caudal,maxY_peak_tType_caudal];
    
    minY_peak_tType_PRC = 0;
    maxY_peak_tType_PRC = p.sizeOfPeak;
    axisBounds_peak_tType_PRC = [0,p.nTrials,minY_peak_tType_PRC,maxY_peak_tType_PRC];
    
    %--------------------------------------------------------------------------
    % LESION
    %--------------------------------------------------------------------------
    
    %**************************************************************************
    % caudal
    %**************************************************************************
    
    % new stimuli
    figs(15) = figure;
    subplot(2,2,1)
    hold on
    plot(peakAct_match_caudal_new_avg(1,:), 'color', 'g')
    plot(peakAct_match_caudal_new_avg(2,:), 'color', 'r')
    plot(peakAct_match_caudal_new_avg(3,:), 'color', 'b')
    
    plot(peakAct_misMatch_caudal_new_avg(1,:), 'g', 'LineStyle', '--')
    plot(peakAct_misMatch_caudal_new_avg(2,:), 'color','r', 'LineStyle', '--')
    plot(peakAct_misMatch_caudal_new_avg(3,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'peakAct, LESION, caudal, NEW'});
    axis(axisBounds_peak_tType_caudal);
    
    subplot(2,2,2); hold on
    plot(totalAct_match_caudal_new_avg(1,:), 'color', 'g')
    plot(totalAct_match_caudal_new_avg(2,:), 'color', 'r')
    plot(totalAct_match_caudal_new_avg(3,:), 'color', 'b')
    
    plot(totalAct_misMatch_caudal_new_avg(1,:), 'g', 'LineStyle', '--')
    plot(totalAct_misMatch_caudal_new_avg(2,:), 'color','r', 'LineStyle', '--')
    plot(totalAct_misMatch_caudal_new_avg(3,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'totalAct, LESION, caudal, NEW'});
    
    % prev stimuli
    subplot(2,2,3); hold on
    plot(peakAct_match_caudal_prev_avg(1,:), 'color', 'g')
    plot(peakAct_match_caudal_prev_avg(2,:), 'color', 'r')
    plot(peakAct_match_caudal_prev_avg(3,:), 'color', 'b')
    
    plot(peakAct_misMatch_caudal_prev_avg(1,:), 'g', 'LineStyle', '--')
    plot(peakAct_misMatch_caudal_prev_avg(2,:), 'color','r', 'LineStyle', '--')
    plot(peakAct_misMatch_caudal_prev_avg(3,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'peakAct, LESION, caudal, PREV'});
    axis(axisBounds_peak_tType_caudal);
    
    subplot(2,2,4); hold on
    plot(totalAct_match_caudal_prev_avg(1,:), 'color', 'g')
    plot(totalAct_match_caudal_prev_avg(2,:), 'color', 'r')
    plot(totalAct_match_caudal_prev_avg(3,:), 'color', 'b')
    
    plot(totalAct_misMatch_caudal_prev_avg(1,:), 'g', 'LineStyle', '--')
    plot(totalAct_misMatch_caudal_prev_avg(2,:), 'color','r', 'LineStyle', '--')
    plot(totalAct_misMatch_caudal_prev_avg(3,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'totalAct, LESION, caudal, PREV'});
    
    saveas(figs(15),[saveFolder, '/pkTA_les_cud'], 'fig');
    saveas(figs(15),[saveFolder, '/pkTA_les_cud'], 'jpg');
    
    
    % no PRC to plot
    
    %--------------------------------------------------------------------------
    % CONTROL
    %--------------------------------------------------------------------------
    
    %**************************************************************************
    % caudal
    %**************************************************************************
    
    % new stimuli
    figs(16) = figure;
    subplot(2,2,1)
    hold on
    plot(peakAct_match_caudal_new_avg(4,:), 'color', 'g')
    plot(peakAct_match_caudal_new_avg(5,:), 'color', 'r')
    plot(peakAct_match_caudal_new_avg(6,:), 'color', 'b')
    
    plot(peakAct_misMatch_caudal_new_avg(4,:), 'g', 'LineStyle', '--')
    plot(peakAct_misMatch_caudal_new_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(peakAct_misMatch_caudal_new_avg(6,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'peakAct, CONTROL, caudal, NEW'});
    axis(axisBounds_peak_tType_caudal);
    
    subplot(2,2,2); hold on
    plot(totalAct_match_caudal_new_avg(4,:), 'color', 'g')
    plot(totalAct_match_caudal_new_avg(5,:), 'color', 'r')
    plot(totalAct_match_caudal_new_avg(6,:), 'color', 'b')
    
    plot(totalAct_misMatch_caudal_new_avg(4,:), 'g', 'LineStyle', '--')
    plot(totalAct_misMatch_caudal_new_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(totalAct_misMatch_caudal_new_avg(6,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'totalAct, CONTROL, caudal, NEW'});
    
    % prev stim
    subplot(2,2,3)
    hold on
    plot(peakAct_match_caudal_prev_avg(4,:), 'color', 'g')
    plot(peakAct_match_caudal_prev_avg(5,:), 'color', 'r')
    plot(peakAct_match_caudal_prev_avg(6,:), 'color', 'b')
    
    plot(peakAct_misMatch_caudal_prev_avg(4,:), 'g', 'LineStyle', '--')
    plot(peakAct_misMatch_caudal_prev_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(peakAct_misMatch_caudal_prev_avg(6,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'peakAct, CONTROL, caudal, PREV'});
    axis(axisBounds_peak_tType_caudal);
    
    subplot(2,2,4); hold on
    plot(totalAct_match_caudal_prev_avg(4,:), 'color', 'g')
    plot(totalAct_match_caudal_prev_avg(5,:), 'color', 'r')
    plot(totalAct_match_caudal_prev_avg(6,:), 'color', 'b')
    
    plot(totalAct_misMatch_caudal_prev_avg(4,:), 'g', 'LineStyle', '--')
    plot(totalAct_misMatch_caudal_prev_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(totalAct_misMatch_caudal_prev_avg(6,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'totalAct, CONTROL, caudal, PREV'});
    
    saveas(figs(16),[saveFolder, '/pkTt_cntr_caud'], 'fig');
    saveas(figs(16),[saveFolder, '/pkTt_cntr_caud'], 'jpg');
    
    %**************************************************************************
    % PRC
    %**************************************************************************
    
    % new stim
    figs(17) = figure;
    subplot(2,2,1)
    hold on
    plot(peakAct_match_PRC_new_avg(4,:), 'color', 'g')
    plot(peakAct_match_PRC_new_avg(5,:), 'color', 'r')
    plot(peakAct_match_PRC_new_avg(6,:), 'color', 'b')
    
    plot(peakAct_misMatch_PRC_new_avg(4,:), 'g', 'LineStyle', '--')
    plot(peakAct_misMatch_PRC_new_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(peakAct_misMatch_PRC_new_avg(6,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'peakAct, PRC, CONTROL, NEW'});
    axis(axisBounds_peak_tType_PRC);
    
    subplot(2,2,2); hold on
    plot(totalAct_match_PRC_new_avg(4,:), 'color', 'g')
    plot(totalAct_match_PRC_new_avg(5,:), 'color', 'r')
    plot(totalAct_match_PRC_new_avg(6,:), 'color', 'b')
    
    plot(totalAct_misMatch_PRC_new_avg(4,:), 'g', 'LineStyle', '--')
    plot(totalAct_misMatch_PRC_new_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(totalAct_misMatch_PRC_new_avg(6,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'totalAct, PRC, CONTROL, NEW'});
    
    % prev stim
    subplot(2,2,3)
    hold on
    plot(peakAct_match_PRC_prev_avg(4,:), 'color', 'g')
    plot(peakAct_match_PRC_prev_avg(5,:), 'color', 'r')
    plot(peakAct_match_PRC_prev_avg(6,:), 'color', 'b')
    
    plot(peakAct_misMatch_PRC_prev_avg(4,:), 'g', 'LineStyle', '--')
    plot(peakAct_misMatch_PRC_prev_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(peakAct_misMatch_PRC_prev_avg(6,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'peakAct, PRC, CONTROL, PREV'});
    axis(axisBounds_peak_tType_PRC);
    
    subplot(2,2,4); hold on
    plot(totalAct_match_PRC_prev_avg(4,:), 'color', 'g')
    plot(totalAct_match_PRC_prev_avg(5,:), 'color', 'r')
    plot(totalAct_match_PRC_prev_avg(6,:), 'color', 'b')
    
    plot(totalAct_misMatch_PRC_prev_avg(4,:), 'g', 'LineStyle', '--')
    plot(totalAct_misMatch_PRC_prev_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(totalAct_misMatch_PRC_prev_avg(6,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('activation');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'totalAct, PRC, CONTROL, PREV'});
    
    
    saveas(figs(17),[saveFolder, '/pkTt_cntr_PRC'], 'fig');
    saveas(figs(17),[saveFolder, '/pkTt_cntr_PRC'], 'jpg');
    
    
    
    %% counting when PRC was used
    
    % previous stim
    figs(18) = figure;
    subplot(1,2,1)
    hold on
    plot(ifPRC_match_prev_avg(4,:), 'color', 'g')
    plot(ifPRC_match_prev_avg(5,:), 'color', 'r')
    plot(ifPRC_match_prev_avg(6,:), 'color', 'b')
    
    plot(ifPRC_misMatch_prev_avg(4,:), 'g', 'LineStyle', '--')
    plot(ifPRC_misMatch_prev_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(ifPRC_misMatch_prev_avg(6,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('how often PRC is used');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'Whether PRC used (5 feats sampled), previous Stim'});
    
    % new stim
    subplot(1,2,2); hold on
    plot(ifPRC_match_new_avg(4,:), 'color', 'g')
    plot(ifPRC_match_new_avg(5,:), 'color', 'r')
    plot(ifPRC_match_new_avg(6,:), 'color', 'r')
    
    plot(ifPRC_misMatch_new_avg(4,:), 'g', 'LineStyle', '--')
    plot(ifPRC_misMatch_new_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(ifPRC_misMatch_new_avg(6,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('how often PRC is used');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'Whether PRC used (5 feats sampled), new Stim'});
    
    saveas(figs(18),[saveFolder, '/ifPRC_newPrev'], 'fig');
    saveas(figs(18),[saveFolder, '/ifPRC_newPrev'], 'jpg');
    
    
    % both stimuli
    figs(19) = figure;
    hold on;
    plot(ifPRC_match_both_avg(4,:), 'color', 'g')
    plot(ifPRC_match_both_avg(5,:), 'color', 'r')
    plot(ifPRC_match_both_avg(6,:), 'color', 'b')
    
    plot(ifPRC_misMatch_both_avg(4,:), 'g', 'LineStyle', '--')
    plot(ifPRC_misMatch_both_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(ifPRC_misMatch_both_avg(6,:), 'color','r', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('how often PRC is used');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'Whether PRC used for trial'});
    
    saveas(figs(19),[saveFolder, '/ifPRC_both'], 'fig');
    saveas(figs(19),[saveFolder, '/ifPRC_both'], 'jpg');
    
    %% plot familiarity differences that were used
    
    %--------------------------------------------------------------------------
    % define axis properties
    %--------------------------------------------------------------------------
    % minY_familDiffs = minY_familDiffs;
    % maxY_familDiffs = maxY_familDiffs;
    axisBounds_familDiffs = [0,p.nTrials,minY_familDiffs,maxY_familDiffs];
    
    
    figs(20) = figure;
    subplot(1,2,1); hold on
    plot(familDiff_used_match_avg(1,:), 'color', 'g')
    plot(familDiff_used_match_avg(2,:), 'color', 'r')
    plot(familDiff_used_match_avg(3,:), 'color', 'b')
    
    plot(familDiff_used_misMatch_avg(1,:), 'g', 'LineStyle', '--')
    plot(familDiff_used_misMatch_avg(2,:), 'color','r', 'LineStyle', '--')
    plot(familDiff_used_misMatch_avg(3,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('familiarity difference');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'familiarity difference during LESION sessions'});
    axis(axisBounds_familDiffs);
    
    subplot(1,2,2); hold on
    plot(familDiff_used_match_avg(4,:), 'color', 'g')
    plot(familDiff_used_match_avg(5,:), 'color', 'r')
    plot(familDiff_used_match_avg(6,:), 'color', 'b')
    
    plot(familDiff_used_misMatch_avg(4,:), 'g', 'LineStyle', '--')
    plot(familDiff_used_misMatch_avg(5,:), 'color','r', 'LineStyle', '--')
    plot(familDiff_used_misMatch_avg(6,:), 'color','b', 'LineStyle', '--')
    
    xlabel('trial');
    ylabel('familiarity difference');
    legend('match,low1', 'match,high', 'match,low2', 'misMatch,low1','misMatch,high','misMatch,low2','Location','best');
    title({'familiarity difference during CONTROL sessions'});
    axis(axisBounds_familDiffs);
    
    saveas(figs(20),[saveFolder, '/familDiff_used'], 'fig');
    saveas(figs(20),[saveFolder, '/familDiff_used'], 'jpg');
    
    %% adjusted dPrime graphs
    
    % %
    % figs(21) = figure;
    % subplot(1,2,1)
    % barweb([dPrime_first_adj_all_mean(2), dPrime_second_adj_all_mean(2) ; dPrime_first_adj_all_mean(1), dPrime_second_adj_all_mean(1)], [dPrime_first_adj_all_err(2), dPrime_second_adj_all_err(2) ; dPrime_first_adj_all_err(1), dPrime_second_adj_all_err(1)], [], {'high', 'low'})
    % xlabel('stim ambiguity');
    % ylabel('d''');
    % legend('first Half', 'Second Half')
    % title({'dPrime of Lesion'})
    %
    % subplot(1,2,2)
    % barweb([dPrime_first_adj_all_mean(4), dPrime_second_adj_all_mean(4) ; dPrime_first_adj_all_mean(3), dPrime_second_adj_all_mean(3)], [dPrime_first_adj_all_err(4), dPrime_second_adj_all_err(4) ; dPrime_first_adj_all_err(3), dPrime_second_adj_all_err(3)], [], {'high', 'low'})
    % xlabel('stim ambiguity');
    % ylabel('d''');
    % legend('first Half', 'Second Half')
    % title({'dPrime of Control'})
    %
    % saveas(figs(21),[saveFolder, '/dPrime_adj'], 'fig');
    % saveas(figs(21),[saveFolder, '/dPrime_adj'], 'jpg');
    
    
    % figs(22) = figure;
    % subplot(1,2,1)
    % barweb([dPrime_first_adj_some_mean(2), dPrime_second_adj_some_mean(2) ; dPrime_first_adj_some_mean(1), dPrime_second_adj_some_mean(1)], [dPrime_first_adj_some_err(2), dPrime_second_adj_some_err(2) ; dPrime_first_adj_some_err(1), dPrime_second_adj_some_err(1)], [], {'high', 'low'})
    % xlabel('stim ambiguity');
    % ylabel('d''');
    % legend('first Half', 'Second Half')
    % title({'dPrime of Lesion'})
    %
    % subplot(1,2,2)
    % barweb([dPrime_first_adj_some_mean(4), dPrime_second_adj_some_mean(4) ; dPrime_first_adj_some_mean(3), dPrime_second_adj_some_mean(3)], [dPrime_first_adj_some_err(4), dPrime_second_adj_some_err(4) ; dPrime_first_adj_some_err(3), dPrime_second_adj_some_err(3)], [], {'high', 'low'})
    % xlabel('stim ambiguity');
    % ylabel('d''');
    % legend('first Half', 'Second Half')
    % title({'dPrime of Control'})
    %
    % saveas(figs(22),[saveFolder, '/dPrime_adj_some'], 'fig');
    % saveas(figs(22),[saveFolder, '/dPrime_adj_some'], 'jpg');
    
    %%
    minY_dPrime = 0;
    maxY_dPrime = 6;
    axisBounds_dPrime = [0,1,minY_dPrime,maxY_dPrime];
    
    
    figs(23) = figure;
    subplot(1,2,1)
    barweb([dP_adj_all_mean(1); dP_adj_all_mean(2); dP_adj_all_mean(3)], [dP_adj_all_sem(1); dP_adj_all_sem(2); dP_adj_all_sem(3)], [], {})
    xlabel('stim ambiguity');
    ylabel('d''');
    legend('LA1', 'HA', 'LA2','Location','southeast')
    title({'dPrime of Lesion'})
    figs(23).CurrentAxes.YLim = [minY_dPrime, maxY_dPrime];
    
    subplot(1,2,2)
    barweb([dP_adj_all_mean(4), dP_adj_all_mean(5), dP_adj_all_mean(6)], [dP_adj_all_sem(4), dP_adj_all_sem(5), dP_adj_all_sem(6)], [], {})
    xlabel('stim ambiguity');
    ylabel('d''');
    legend('LA1', 'HA', 'LA2','Location','southeast')
    title({'dPrime of Control'})
    figs(23).CurrentAxes.YLim = [minY_dPrime, maxY_dPrime];
    
    saveas(figs(23),[saveFolder, '/dPrime_adj'], 'fig');
    saveas(figs(23),[saveFolder, '/dPrime_adj'], 'jpg');
    
    
    % d' of comparison trials only
    figs(24) = figure;
    subplot(1,2,1)
    barweb([dP_all_comp_mean(1); dP_all_comp_mean(2); dP_all_comp_mean(3)], [dP_all_comp_sem(1); dP_all_comp_sem(2); dP_all_comp_sem(3)], [], {})
    xlabel('stim ambiguity');
    ylabel('d''');
    legend('LA1', 'HA', 'LA2','Location','southeast')
    title({'dPrime of Lesion, comp trials only'})
    figs(24).CurrentAxes.YLim = [minY_dPrime, maxY_dPrime];
    
    
    subplot(1,2,2)
    barweb([dP_all_comp_mean(4), dP_all_comp_mean(5), dP_all_comp_mean(6)], [dP_all_comp_sem(4), dP_all_comp_sem(5), dP_all_comp_sem(6)], [], {})
    xlabel('stim ambiguity');
    ylabel('d''');
    legend('LA1', 'HA', 'LA2','Location','southeast')
    title({'dPrime of Control, comp trials only'})
    figs(24).CurrentAxes.YLim = [minY_dPrime, maxY_dPrime];
    
    
    saveas(figs(24),[saveFolder, '/dPrime_adj_comp'], 'fig');
    saveas(figs(24),[saveFolder, '/dPrime_adj_comp'], 'jpg');
    
    %%
else
    barvalues = [dP_all_comp_mean(4), dP_all_comp_mean(5), dP_all_comp_mean(6)];
    errors = [dP_all_comp_sem(4), dP_all_comp_sem(5), dP_all_comp_sem(6)];
    
    figs(25) = figure;
    handles = barweb(barvalues...
        , errors...
        , [] ...
        , {'Low 1              High              Low 2'}...
        , {'Simulation 2'}...
        , {'Interference'}...
        , {'d'''}...
        , [rgb('Moccasin') ; rgb('Moccasin') ; rgb('Moccasin')]);
    set(gca,'fontsize',30)
    figs(25).CurrentAxes.YLim = [0, 6];
    set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
    
    x1 = handles.bars(1).XOffset;
    x2 = handles.bars(2).XOffset;
    x3 = handles.bars(3).XOffset;
    
    barvalues = [dP_all_comp_mean(1); dP_all_comp_mean(2); dP_all_comp_mean(3)];
    errors = [dP_all_comp_sem(1); dP_all_comp_sem(2); dP_all_comp_sem(3)];
    
    hold on
    hand1 = errorbar([1+x1, 1+x2, 1+x3], barvalues, errors, '-kx', 'MarkerSize', 10,'linewidth', 2);
    
    %     legend('Low Ambiguity 1', 'High Ambiguity', 'Low Ambiguity 2', 'Lesion', 'Location','NorthWest');
    legend('Control', 'Location','NorthWest');
    legend BOXOFF
    
    saveas(figs(25),[saveFolder, '/dPrime_paper'], 'fig');
    saveas(figs(25),[saveFolder, '/dPrime_paper'], 'jpg');
    
end

end

