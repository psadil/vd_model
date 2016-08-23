function dPrimePredictions = calcDPrime(firstRat, lastRat, folderName)

% analyze familiarity differences
% analyze familiarity differences

% now called directly from create_sim

% currently, taking the maximum familiarity differences at each trial
% produces the set of most desirable results
% folderName must be a string (without / on either side)

% scrsz = get(groot, 'ScreenSize');

saveFolder = [pwd,'/graphsAndSession/', folderName];

% load sample data file to get size of parameters
fileName = [saveFolder, '/Session', num2str(1), '_Rat', num2str(1)];
load(fileName)

numRats = lastRat-firstRat+1;

tType = zeros(numRats,4,p.nTrials);

familDifferences_misMatch_caudal = zeros(numRats,4,p.nTrials);
familDifferences_match_caudal = zeros(numRats,4,p.nTrials);
familDifferences_match_PRC = zeros(numRats,4,p.nTrials);
familDifferences_misMatch_PRC = zeros(numRats,4,p.nTrials);

% familDiff_caudal = zeros(numRats,4,p.nTrials);
% familDiff_PRC = zeros(numRats,4,p.nTrials);

answer = zeros(numRats,4,p.nTrials);
correct = zeros(numRats,4,p.nTrials);
acc_firstHalf = zeros(numRats,4);
acc_secondHalf = zeros(numRats,4);
acc_match_first = zeros(numRats,4);
acc_misMatch_first = zeros(numRats,4);
acc_match_second = zeros(numRats,4);
acc_misMatch_second = zeros(numRats,4);

dPrime_first = zeros(numRats,4);
dPrime_second = zeros(numRats,4);

hitRate_first = zeros(numRats,4);
hitRate_second = zeros(numRats,4);
FARate_first = zeros(numRats,4);
FARate_second = zeros(numRats,4);

meanSelectivity_caudal_prev = zeros(numRats,4,p.nTrials);
meanSelectivity_caudal_new = zeros(numRats,4,p.nTrials);
meanSelectivity_PRC_prev = zeros(numRats,4,p.nTrials);
meanSelectivity_PRC_new = zeros(numRats,4,p.nTrials);

familDiff_used = zeros(numRats,4,p.nTrials);

meanSelectivity_caudal_prev_misMatch = zeros(numRats,4,p.nTrials);
meanSelectivity_caudal_prev_match = zeros(numRats,4,p.nTrials);
meanSelectivity_caudal_new_misMatch = zeros(numRats,4,p.nTrials);
meanSelectivity_caudal_new_match = zeros(numRats,4,p.nTrials);

meanSelectivity_PRC_prev_misMatch = zeros(numRats,4,p.nTrials);
meanSelectivity_PRC_prev_match = zeros(numRats,4,p.nTrials);
meanSelectivity_PRC_new_misMatch = zeros(numRats,4,p.nTrials);
meanSelectivity_PRC_new_match = zeros(numRats,4,p.nTrials);

comparedFeat = zeros(numRats,4,p.nTrials,p.numGrids_Caudal);

trial_misMatch_caudal = zeros(4,p.nTrials);
trial_match_caudal = zeros(4,p.nTrials);
trial_misMatch_PRC = zeros(4,p.nTrials);
trial_match_PRC = zeros(4,p.nTrials);

trial_misMatch_caudal_if = zeros(numRats,4,p.nTrials);
trial_match_caudal_if = zeros(numRats,4,p.nTrials);
trial_misMatch_PRC_if = zeros(numRats,4,p.nTrials);
trial_match_PRC_if = zeros(numRats,4,p.nTrials);

gridUsed_one_total = zeros(4,p.nTrials);
gridUsed_one = zeros(numRats,4,p.nTrials);

gridUsed_one_misMatch_total = zeros(4,p.nTrials);
gridUsed_one_misMatch = zeros(numRats,4,p.nTrials);

gridUsed_one_match_total = zeros(4,p.nTrials);
gridUsed_one_match = zeros(numRats,4,p.nTrials);

selectivity_caudal_prev_match_grid1 = zeros(numRats,4,p.nTrials);
selectivity_caudal_prev_misMatch_grid1 = zeros(numRats,4,p.nTrials);
selectivity_caudal_new_match_grid1 = zeros(numRats,4,p.nTrials);
selectivity_caudal_new_misMatch_grid1 = zeros(numRats,4,p.nTrials);


fixationByComparison_prev = zeros(numRats,4,p.nTrials);
fixationByComparison_new = zeros(numRats,4,p.nTrials);

fixByComp_prev_mean = zeros(numRats,4,p.nTrials,max(p.maxFixations));
fixByComp_new_mean = zeros(numRats,4,p.nTrials,max(p.maxFixations));


peakAct = zeros(numRats,4,p.nTrials,2);
totalAct = zeros(numRats,4,p.nTrials,2);

featSamplePerTrial = zeros(numRats,4,p.nTrials);
featSamplePerTrial_misMatch = zeros(numRats,4,p.nTrials);
featSamplePerTrial_match = zeros(numRats,4,p.nTrials);
featSamplePerTrial_misMatch_count = zeros(4,p.nTrials);
featSamplePerTrial_match_count = zeros(4,p.nTrials);

ifPRC_both = zeros(numRats,4,2,p.nTrials);

act_peak_prev_init = zeros(numRats,4,2,p.nTrials);
act_peak_prev_fin = zeros(numRats,4,2,p.nTrials);
act_peak_new_init = zeros(numRats,4,2,p.nTrials);
act_peak_new_fin = zeros(numRats,4,2,p.nTrials);
act_total_prev_init = zeros(numRats,4,2,p.nTrials);
act_total_prev_fin = zeros(numRats,4,2,p.nTrials);
act_total_new_init = zeros(numRats,4,2,p.nTrials);
act_total_new_fin = zeros(numRats,4,2,p.nTrials);

fixations_total = zeros(numRats,4,p.nTrials);

for rat = firstRat:lastRat
    for session = 1:4
        
        fileName = [saveFolder, '/Session', num2str(session), '_Rat', num2str(rat)];
        load(fileName)
        fprintf ('\nloading rat %d, session %d', rat, session);
        
        tType(rat,session,:) = p.tType;
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
        if session > 2
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
        
        if session > 2
            
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
        
        
        
        % not sure anymore what this does, something about arranging number
        % of fixations...
        if any(session == [1,3])
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
            
            if p.tType(trial) == 1
                
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
                
                
                
                
                if session > 2 && (all(p.usePRC(:,trial)))
                    %                     meanSelectivity_PRC_prev_misMatch(rat,session,trial) = p.meanSelectivity_PRC_prev(trial);
                    %                     meanSelectivity_PRC_new_misMatch(rat,session,trial) = p.meanSelectivity_PRC_new(trial);
                    
                    familDifferences_misMatch_PRC(rat,session,trial) = p.familDiff_PRC(trial);
                    trial_misMatch_PRC(session,trial) = trial_misMatch_PRC(session,trial) + 1;
                    trial_misMatch_PRC_if(rat,session,trial) = 1;
                end
                
                %--------------------------------------------------------------
                % match Trials
                %--------------------------------------------------------------
                
            elseif p.tType(trial) == 2
                
                
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
                if session > 2 && (all(p.usePRC(:,trial)))
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
        
        
        % adjust by adding to all
        hitRate_first_adj_all(rat,session) = (sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==1)) + 1/(2*p.nTrials))/(sum(p.tType(1:p.nTrials/2)==1)+1);  % a 'yes' (mismatch judgement) on trials that were mismatches (ie, p.tType==1)
        FARate_first_adj_all(rat,session) = (sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==2)) + 1/(2*p.nTrials))/(sum(p.tType(1:p.nTrials/2)==2)+1); % a 'yes' on matching trials
        
        hitRate_second_adj_all(rat,session) = (sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==1)) + 1/(2*p.nTrials))/(sum(p.tType(p.nTrials/2+1:end)==1)+1);  % a 'yes' (mismatch judgement) on trials that were mismatches (ie, p.tType==1)
        FARate_second_adj_all(rat,session) = (sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==2)) + 1/(2*p.nTrials))/(sum(p.tType(p.nTrials/2+1:end)==2)+1); % a 'yes' on matching trials
        
        % adjust by adding to only trials
        hitRate_first_adj_some(rat,session) = (sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==1)) + 1/(2*p.nTrials))/(sum(p.tType(1:p.nTrials/2)==1)+1);  % a 'yes' (mismatch judgement) on trials that were mismatches (ie, p.tType==1)
        if hitRate_first_adj_some(rat,session) == 1
            hitRate_first_adj_some(rat,session) = hitRate_first_adj_some(rat,session) - 1/(2*p.nTrials);
        elseif hitRate_first_adj_some(rat,session) == 0
            hitRate_first_adj_some(rat,session) = hitRate_first_adj_some(rat,session) + 1/(2*p.nTrials);
        end
        
        FARate_first_adj_some(rat,session) = (sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==2)) + 1/(2*p.nTrials))/(sum(p.tType(1:p.nTrials/2)==2)+1); % a 'yes' on matching trials
        if FARate_first_adj_some(rat,session) == 1
            FARate_first_adj_some(rat,session) = FARate_first_adj_some(rat,session) - 1/(2*p.nTrials);
        elseif FARate_first_adj_some(rat,session) == 0
            FARate_first_adj_some(rat,session) = FARate_first_adj_some(rat,session) + 1/(2*p.nTrials);
        end
        
        hitRate_second_adj_some(rat,session) = (sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==1)) + 1/(2*p.nTrials))/(sum(p.tType(p.nTrials/2+1:end)==1)+1);  % a 'yes' (mismatch judgement) on trials that were mismatches (ie, p.tType==1)
        if hitRate_second_adj_some(rat,session) == 1
            hitRate_second_adj_some(rat,session) = hitRate_second_adj_some(rat,session) - 1/(2*p.nTrials);
        elseif hitRate_second_adj_some(rat,session) == 0
            hitRate_second_adj_some(rat,session) = hitRate_second_adj_some(rat,session) + 1/(2*p.nTrials);
        end
        
        FARate_second_adj_some(rat,session) = (sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==2)) + 1/(2*p.nTrials))/(sum(p.tType(p.nTrials/2+1:end)==2)+1); % a 'yes' on matching trials
        if FARate_second_adj_some(rat,session) == 1
            FARate_second_adj_some(rat,session) = FARate_second_adj_some(rat,session) - 1/(2*p.nTrials);
        elseif FARate_second_adj_some(rat,session) == 0
            FARate_second_adj_some(rat,session) = FARate_second_adj_some(rat,session) + 1/(2*p.nTrials);
        end
        
        
        acc_firstHalf(rat,session) = p.Acc_firstHalf;
        acc_secondHalf(rat,session) = p.Acc_secondHalf;
        
        acc_match_first(rat,session) = sum(answer(rat,session,find(tType(rat,session,(1:36))==2))==0)/(sum(tType(rat,session,1:36)==2));
        acc_match_second(rat,session) = sum(answer(rat,session,find(tType(rat,session,(36:p.nTrials))==2))==0)/(sum(tType(rat,session,37:p.nTrials)==2));
        acc_misMatch_first(rat,session) = sum(answer(rat,session,find(tType(rat,session,(1:36))==1))==1)/(sum(tType(rat,session,1:36)==1));
        acc_misMatch_second(rat,session) = sum(answer(rat,session,find(tType(rat,session,(36:p.nTrials))==1))==1)/(sum(tType(rat,session,37:p.nTrials)==1));
        
        
        %         % divide selectivity by peak activation
        %         peakAct(rat,session,:,:) = p.peak_act;
        %         totalAct(rat,session,:,:) = p.totalAct;
    end
end


%% looking at more raw hit and FA rats

hitRate_first_mean = mean(hitRate_first);
hitRate_second_mean = mean(hitRate_second);
FARate_first_mean = mean(FARate_first);
FARate_second_mean = mean(FARate_second);

hitRate_first_std = std(hitRate_first);
hitRate_first_sem = hitRate_first_std/sqrt(numRats);

FARate_first_std = std(FARate_first);
FARate_first_sem = FARate_first_std/sqrt(numRats);

hitRate_second_std = std(hitRate_second);
hitRate_second_sem = hitRate_second_std/sqrt(numRats);

FARate_second_std = std(FARate_second);
FARate_second_sem = FARate_second_std/sqrt(numRats);

dPrime_first_raw = norminv(hitRate_first_mean) - norminv(FARate_first_mean);
dPrime_second_raw = norminv(hitRate_second_mean) - norminv(FARate_second_mean);

dPrime_first_err = abs(norminv(hitRate_first_sem) - norminv(FARate_first_sem));
dPrime_second_err = abs(norminv(hitRate_second_sem) - norminv(FARate_second_sem));


% work out adjusted, adding to all

hitRate_first_adj_all_mean = mean(hitRate_first_adj_all);
hitRate_second_adj_all_mean = mean(hitRate_second_adj_all);

FARate_first_adj_all_mean = mean(FARate_first_adj_all);
FARate_second_adj_all_mean = mean(FARate_second_adj_all);

hitRate_first_adj_all_std = std(hitRate_first_adj_all);
hitRate_first_adj_all_sem = hitRate_first_adj_all_std/sqrt(numRats);

FARate_first_adj_all_std = std(FARate_first_adj_all);
FARate_first_adj_all_sem = FARate_first_adj_all_std/sqrt(numRats);

hitRate_second_adj_all_std = std(hitRate_second_adj_all);
hitRate_second_adj_all_sem = hitRate_second_adj_all_std/sqrt(numRats);

FARate_second_adj_all_std = std(FARate_second_adj_all);
FARate_second_adj_all_sem = FARate_second_adj_all_std/sqrt(numRats);

dPrime_first_adj_all = norminv(hitRate_first_adj_all_mean) - norminv(FARate_first_adj_all_mean);
dPrime_second_adj_all = norminv(hitRate_second_adj_all_mean) - norminv(FARate_second_adj_all_mean);

dPrime_first_adj_all_err = abs(norminv(hitRate_first_adj_all_sem) - norminv(FARate_first_adj_all_sem));
dPrime_second_adj_all_err = abs(norminv(hitRate_second_adj_all_sem) - norminv(FARate_second_adj_all_sem));

% adjusting only rates equal to 0 or 1
hitRate_first_adj_some_mean = mean(hitRate_first_adj_some);
hitRate_second_adj_some_mean = mean(hitRate_second_adj_some);

FARate_first_adj_some_mean = mean(FARate_first_adj_some);
FARate_second_adj_some_mean = mean(FARate_second_adj_some);

hitRate_first_adj_some_std = std(hitRate_first_adj_some);
hitRate_first_adj_some_sem = hitRate_first_adj_some_std/sqrt(numRats);

FARate_first_adj_some_std = std(FARate_first_adj_some);
FARate_first_adj_some_sem = FARate_first_adj_some_std/sqrt(numRats);

hitRate_second_adj_some_std = std(hitRate_second_adj_some);
hitRate_second_adj_some_sem = hitRate_second_adj_some_std/sqrt(numRats);

FARate_second_adj_some_std = std(FARate_second_adj_some);
FARate_second_adj_some_sem = FARate_second_adj_some_std/sqrt(numRats);


dPrime_first_adj_some = norminv(hitRate_first_adj_some_mean) - norminv(FARate_first_adj_some_mean);
dPrime_second_adj_some = norminv(hitRate_second_adj_some_mean) - norminv(FARate_second_adj_some_mean);


dPrime_first_adj_some_err = abs(norminv(hitRate_first_adj_some_sem) - norminv(FARate_first_adj_some_sem));
dPrime_second_adj_some_err = abs(norminv(hitRate_second_adj_some_sem) - norminv(FARate_second_adj_some_sem));



%% Accuracy
acc_match_first_avg = squeeze(mean(acc_match_first,1));
acc_misMatch_first_avg = squeeze(mean(acc_misMatch_first,1));
acc_match_second_avg = squeeze(mean(acc_match_second,1));
acc_misMatch_second_avg = squeeze(mean(acc_misMatch_second,1));

acc_match_bar = [acc_match_first_avg; acc_match_second_avg];
acc_misMatch_bar = [acc_misMatch_first_avg; acc_misMatch_second_avg];

dPrime_first_avg = mean(~isinf(dPrime_first),1);
dPrime_second_avg = mean(~isinf(dPrime_second),1);
acc_firstHalf_avg = mean(acc_firstHalf,1);
acc_secondHalf_avg = mean(acc_secondHalf,1);



%% Now, actually calculate dPrime as needed
% dPrime_first_raw and dPrime_second_raw => (caudal_LA caudal_HA, PRC_LA, PRC_HA)

dPrimeTallied(1,1:4) = dPrime_first_adj_all;
dPrimeTallied(2,1:4) = dPrime_second_adj_all;

dPrimePredictions = dPrimeTallied(2,:) - dPrimeTallied(1,:);

end

