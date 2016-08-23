% old plotFamilDiff routines

% broken down by grid
if p.comparedFeat(trial,1)
    selectivity_caudal_prev_misMatch_grid1(rat,session,trial) = p.meanSelectivity_caudal_prev(trial,1);
    selectivity_caudal_new_misMatch_grid1(rat,session,trial) = p.meanSelectivity_caudal_new(trial,1);
    
    gridUsed_one_misMatch(rat,session,trial) = 1;
    gridUsed_one_misMatch_total(session,trial) = gridUsed_one_misMatch_total(session,trial) + 1;
elseif p.comparedFeat(trial,2)
    selectivity_comparedFive_caudal_prev_misMatch2(rat,session,trial) = p.meanSelectivity_caudal_prev(trial,2);
elseif p.comparedFeat(trial,3)
    selectivity_comparedFive_caudal_prev_misMatch3(rat,session,trial) = p.meanSelectivity_caudal_prev(trial,3);
elseif p.comparedFeat(trial,4)
    selectivity_comparedFive_caudal_prev_misMatch4(rat,session,trial) = p.meanSelectivity_caudal_prev(trial,4);
elseif p.comparedFeat(trial,5)
    selectivity_comparedFive_caudal_prev_misMatch5(rat,session,trial) = p.meanSelectivity_caudal_prev(trial,5);
end


%%

%% Selectivity Absolute of Grid 1

%--------------------------------------------------------------------------
% Caudal
%--------------------------------------------------------------------------

grid1_misMatch_temp = squeeze(sum(selectivity_caudal_prev_misMatch_grid1,1));
grid1_misMatch_mean = grid1_misMatch_temp(:,:) ./ gridUsed_one_misMatch_total(:,:);

grid1_match_temp = squeeze(sum(selectivity_caudal_prev_match_grid1,1));
grid1_match_mean = grid1_match_temp(:,:) ./ gridUsed_one_match_total(:,:);

grid1_misMatch_temp_new = squeeze(sum(selectivity_caudal_new_misMatch_grid1,1));
grid1_misMatch_mean_new = grid1_misMatch_temp_new(:,:) ./ gridUsed_one_misMatch_total(:,:);

grid1_match_temp_new = squeeze(sum(selectivity_caudal_new_match_grid1,1));
grid1_match_mean_new = grid1_match_temp_new(:,:) ./ gridUsed_one_match_total(:,:);

%--------------------------------------------------------------------------
% PRC
%--------------------------------------------------------------------------

selec_PRC_prev_match_temp = squeeze(sum(meanSelectivity_PRC_prev_match(:,:,:),1));
selec_PRC_prev_match = selec_PRC_prev_match_temp ./ trial_match_PRC(:,:);

selec_PRC_prev_misMatch_temp = squeeze(sum(meanSelectivity_PRC_prev_misMatch(:,:,:),1));
selec_PRC_prev_misMatch = selec_PRC_prev_misMatch_temp ./ trial_misMatch_PRC(:,:);

selec_PRC_new_match_temp = squeeze(sum(meanSelectivity_PRC_new_match(:,:,:),1));
selec_PRC_new_match = selec_PRC_new_match_temp ./ trial_match_PRC(:,:);

selec_PRC_new_misMatch_temp = squeeze(sum(meanSelectivity_PRC_new_misMatch(:,:,:),1));
selec_PRC_new_misMatch = selec_PRC_new_misMatch_temp ./ trial_misMatch_PRC(:,:);


%%
% %% for stats
%
% %--------------------------------------------------------------------------
% % caudal
% familDifferences_match_caudal_stat_first_temp = familDifferences_match_caudal(:,1:2,1:36);
% familDifferences_match_caudal_low_stat_first = familDifferences_match_caudal_stat_first_temp(trial_match_caudal_if(:,1,1:36)==1);
% familDifferences_match_caudal_high_stat_first = familDifferences_match_caudal_stat_first_temp(trial_match_caudal_if(:,2,1:36)==1);
%
% familDifferences_match_caudal_stat_second_temp = familDifferences_match_caudal(:,1:2,37:end);
% familDifferences_match_caudal_low_stat_second = familDifferences_match_caudal_stat_second_temp(trial_match_caudal_if(:,1,37:end)==1);
% familDifferences_match_caudal_high_stat_second = familDifferences_match_caudal_stat_second_temp(trial_match_caudal_if(:,2,37:end)==1);
%
% familDifferences_misMatch_caudal_stat_first_temp = familDifferences_misMatch_caudal(:,1:2,1:36);
% familDifferences_misMatch_caudal_low_stat_first = familDifferences_misMatch_caudal_stat_first_temp(trial_misMatch_caudal_if(:,1,1:36)==1);
% familDifferences_misMatch_caudal_high_stat_first = familDifferences_misMatch_caudal_stat_first_temp(trial_misMatch_caudal_if(:,2,1:36)==1);
%
% familDifferences_misMatch_caudal_stat_second_temp = familDifferences_misMatch_caudal(:,1:2,37:end);
% familDifferences_misMatch_caudal_low_stat_second = familDifferences_misMatch_caudal_stat_second_temp(trial_misMatch_caudal_if(:,1,37:end)==1);
% familDifferences_misMatch_caudal_high_stat_second = familDifferences_misMatch_caudal_stat_second_temp(trial_misMatch_caudal_if(:,2,37:end)==1);
%
% %--------------------------------------------------------------------------
% % PRC
% familDifferences_match_PRC_stat_first_temp = familDifferences_match_PRC(:,3:4,1:36);
% familDifferences_match_PRC_low_stat_first = familDifferences_match_PRC_stat_first_temp(trial_match_PRC_if(:,3,1:36)==1);
% familDifferences_match_PRC_high_stat_first = familDifferences_match_PRC_stat_first_temp(trial_match_PRC_if(:,4,1:36)==1);
%
% familDifferences_match_PRC_stat_second_temp = familDifferences_match_PRC(:,3:4,37:end);
% familDifferences_match_PRC_low_stat_second = familDifferences_match_PRC_stat_second_temp(trial_match_PRC_if(:,3,37:end)==1);
% familDifferences_match_PRC_high_stat_second = familDifferences_match_PRC_stat_second_temp(trial_match_PRC_if(:,4,37:end)==1);
%
% familDifferences_misMatch_PRC_stat_first_temp = familDifferences_misMatch_PRC(:,3:4,1:36);
% familDifferences_misMatch_PRC_low_stat_first = familDifferences_misMatch_PRC_stat_first_temp(trial_misMatch_PRC_if(:,3,1:36)==1);
% familDifferences_misMatch_PRC_high_stat_first = familDifferences_misMatch_PRC_stat_first_temp(trial_misMatch_PRC_if(:,4,1:36)==1);
%
% familDifferences_misMatch_PRC_stat_second_temp = familDifferences_misMatch_PRC(:,3:4,37:end);
% familDifferences_misMatch_PRC_low_stat_second = familDifferences_misMatch_PRC_stat_second_temp(trial_misMatch_PRC_if(:,3,37:end)==1);
% familDifferences_misMatch_PRC_high_stat_second = familDifferences_misMatch_PRC_stat_second_temp(trial_misMatch_PRC_if(:,4,37:end)==1);
%
%
% % outPut = [familDifferences_match_caudal_low_stat_first, familDifferences_match_caudal_low_stat_second', familDifferences_match_caudal_low_stat_first', familDifferences_match_caudal_low_stat_first', familDifferences_match_caudal_low_stat_first', familDifferences_match_caudal_low_stat_first']


%%

%% aboslute famil

%--------------------------------------------------------------------------
% LESION
%--------------------------------------------------------------------------

% caudal, prev
figs(4) = figure;
subplot(1,2,1)
hold on
plot(grid1_match_mean(1,:), 'color', 'g')
plot(grid1_match_mean(2,:), 'color', 'r')
plot(grid1_misMatch_mean(1,:), 'g', 'LineStyle', '--')
plot(grid1_misMatch_mean(2,:), 'color','r', 'LineStyle', '--')
legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
title({'Absolute Familiarity, caudal, grid1, LESION'});

% caudal, new
subplot(1,2,2); hold on
plot(grid1_match_mean_new(1,:), 'color', 'g')
plot(grid1_match_mean_new(2,:), 'color', 'r')
plot(grid1_misMatch_mean_new(1,:), 'g', 'LineStyle', '--')
plot(grid1_misMatch_mean_new(2,:), 'color','r', 'LineStyle', '--')
legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
title({'Absolute Familiarity,NEW, caudal, grid1, LESION'});

saveas(figs(4),[saveFolder, '/absolFamil_caud_les'], 'fig');
saveas(figs(4),[saveFolder, '/absolFamil_caud_les'], 'jpg');

% no PRC to plot

%--------------------------------------------------------------------------
% next, plot CONTROL sessions
%--------------------------------------------------------------------------

% caudal, prev
figs(5) = figure;
subplot(1,2,1)
hold on
plot(grid1_match_mean(3,:), 'color', 'g')
plot(grid1_match_mean(4,:), 'color', 'r')
plot(grid1_misMatch_mean(3,:), 'g', 'LineStyle', '--')
plot(grid1_misMatch_mean(4,:), 'color','r', 'LineStyle', '--')
legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
title({'Absolute Familiarity, caudal, grid1, CONTROl'});

% caudal, new
subplot(1,2,2); hold on
plot(grid1_match_mean_new(3,:), 'color', 'g')
plot(grid1_match_mean_new(4,:), 'color', 'r')
plot(grid1_misMatch_mean_new(3,:), 'g', 'LineStyle', '--')
plot(grid1_misMatch_mean_new(4,:), 'color','r', 'LineStyle', '--')
legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
title({'Absolute Familiarity,NEW, caudal, grid1, CONTROL'});

saveas(figs(5),[saveFolder, '/absolFamil_caud_contr'], 'fig');
saveas(figs(5),[saveFolder, '/absolFamil_caud_contr'], 'jpg');

% PRC
figs(6) = figure;
subplot(1,2,1)
hold on
plot(selec_PRC_prev_match(3,:), 'color', 'g')
plot(selec_PRC_prev_match(4,:), 'color', 'r')
plot(selec_PRC_prev_misMatch(3,:), 'g', 'LineStyle', '--')
plot(selec_PRC_prev_misMatch(4,:), 'color','r', 'LineStyle', '--')
legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
title({'Absolute Familiarity, PRC'});

subplot(1,2,2); hold on
plot(selec_PRC_new_match(3,:), 'color', 'g')
plot(selec_PRC_new_match(4,:), 'color', 'r')
plot(selec_PRC_new_misMatch(3,:), 'g', 'LineStyle', '--')
plot(selec_PRC_new_misMatch(4,:), 'color','r', 'LineStyle', '--')
legend('low,match', 'high,match', 'low,misMatch', 'high,misMatch');
title({'Absolute Familiarity,NEW, PRC'});

saveas(figs(6),[saveFolder, '/absolFamil_PRC'], 'fig');
saveas(figs(6),[saveFolder, '/absolFamil_PRC'], 'jpg');

%%

% broken down by grid
                if p.comparedFeat(trial,1)
                    selectivity_caudal_prev_match_grid1(rat,session,trial) = p.meanSelectivity_caudal_prev(trial,1);
                    selectivity_caudal_new_match_grid1(rat,session,trial) = p.meanSelectivity_caudal_new(trial,1);
                    
                    gridUsed_one_match(rat,session,trial) = 1;
                    gridUsed_one_match_total(session,trial) = gridUsed_one_match_total(session,trial) + 1;
                    
                elseif p.comparedFeat(trial,2)
                    selectivity_comparedFive_caudal_prev_match2(rat,session,trial) = p.meanSelectivity_caudal_prev(trial,2);
                elseif p.comparedFeat(trial,3)
                    selectivity_comparedFive_caudal_prev_match3(rat,session,trial) = p.meanSelectivity_caudal_prev(trial,3);
                elseif p.comparedFeat(trial,4)
                    selectivity_comparedFive_caudal_prev_match4(rat,session,trial) = p.meanSelectivity_caudal_prev(trial,4);
                elseif p.comparedFeat(trial,5)
                    selectivity_comparedFive_caudal_prev_match5(rat,session,trial) = p.meanSelectivity_caudal_prev(trial,5);
                end

                
                %%
                dPrime_first(rat,session) = norminv(hitRate_first(rat,session))-norminv(FARate_first(rat,session));
        dPrime_second(rat,session) = norminv(hitRate_second(rat,session))-norminv(FARate_second(rat,session));
                
                
                dPrime_first_avg = mean(~isinf(dPrime_first),1);
dPrime_second_avg = mean(~isinf(dPrime_second),1);
acc_firstHalf_avg = mean(acc_firstHalf,1);
acc_secondHalf_avg = mean(acc_secondHalf,1);

% dPrime_first_match = dPrime_first_avg.*(mean(tTypeProp_match(:,1:36),2))';
% dPrime_second_match = dPrime_second_avg.*(mean(tTypeProp_match(:,37:end),2))';
% dPrime_first_misMatch = dPrime_first_avg.*(mean(tTypeProp_misMatch(:,1:36),2))';
% dPrime_second_misMatch = dPrime_second_avg.*(mean(tTypeProp_misMatch(:,37:end),2))';

dPrime_caudal = [dPrime_first_avg(1:2);dPrime_second_avg(1:2)];
dPrime_PRC = [dPrime_first_avg(3:4);dPrime_second_avg(3:4)];


dPrime_bar = [dPrime_caudal, dPrime_PRC];
dPrime_firstSecond = [dPrime_first_avg; dPrime_second_avg];

%%

% %% peak and total activations
% 
% % caudal
% figs(10) = figure;
% subplot(1,2,1)
% hold on
% plot(peakAct_mean(1,:,1), 'color', 'g')
% plot(peakAct_mean(2,:,1), 'color', 'r')
% plot(peakAct_mean(3,:,1), 'g', 'LineStyle', '--')
% plot(peakAct_mean(4,:,1), 'color','r', 'LineStyle', '--')
% xlabel('trial');
% ylabel('activation');
% legend('lesion,low', 'lesion,high', 'control,low', 'control,high');
% title({'peakAct, caudal, NEW'});
% 
% subplot(1,2,2); hold on
% plot(totalAct_mean(1,:,1), 'color', 'g')
% plot(totalAct_mean(2,:,1), 'color', 'r')
% plot(totalAct_mean(3,:,1), 'g', 'LineStyle', '--')
% plot(totalAct_mean(4,:,1), 'color','r', 'LineStyle', '--')
% xlabel('trial');
% ylabel('activation');
% legend('lesion,low', 'lesion,high', 'control,low', 'control,high');
% title({'totalAct, caudal, NEW'});
% 
% saveas(figs(10),[saveFolder, '/peakOtherActs_caudal'], 'fig');
% saveas(figs(10),[saveFolder, '/peakOtherActs_caudal'], 'jpg');
% 
% % PRC
% figs(11) = figure;
% subplot(1,2,1)
% hold on
% plot(peakAct_mean(3,:,2), 'color', 'g')
% plot(peakAct_mean(4,:,2), 'color', 'r')
% xlabel('trial');
% ylabel('activation');
% legend('control,low', 'control,high');
% title({'peakAct, PRC, NEW'});
% 
% subplot(1,2,2); hold on
% plot(totalAct_mean(3,:,2), 'g', 'LineStyle', '--')
% plot(totalAct_mean(4,:,2), 'color','r', 'LineStyle', '--')
% xlabel('trial');
% ylabel('activation');
% legend('control,low', 'control,high');
% title({'totalAct, PRC, NEW'});
% 
% saveas(figs(11),[saveFolder, '/peakOtherActs_PRC'], 'fig');
% saveas(figs(11),[saveFolder, '/peakOtherActs_PRC'], 'jpg');
