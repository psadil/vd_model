function dPrimePredictions = calcDPrime(firstRat, lastRat, folderName)

% analyze familiarity differences

% now called directly from create_sim

% currently, taking the maximum familiarity differences at each trial
% produces the set of most desirable results
% folderName must be a string (without / on either side)

% folderName = '1encod_p6Ratio1p2_1sampVar5_0noise_100train_20Max25_5peak_100rows_A01_B4_forceSamp_NOshrinkingLearn_etaExp5p_Gexpp8_stim4Diff_NOfives';

% scrsz = get(groot, 'ScreenSize');

saveFolder = [pwd,'/graphsAndSession/', folderName];

% load sample data file to get size of parameters
fileName = [saveFolder, '/Session', num2str(1), '_Rat', num2str(1)];
load(fileName)

numRats = lastRat-firstRat+1;

hitRate_first = zeros(numRats,4);
hitRate_second = zeros(numRats,4);
FARate_first = zeros(numRats,4);
FARate_second = zeros(numRats,4);

for rat = firstRat:lastRat
    for session = 1:4
        
        fileName = [saveFolder, '/Session', num2str(session), '_Rat', num2str(rat)];
        load(fileName)
        
        
        %------------------------------------------------------------------
        % tally results
        %------------------------------------------------------------------
        
        yes1 = sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==1));
        yes1fa = sum(p.answer(1:p.nTrials/2).*(p.tType(1:p.nTrials/2)==2));
        yes2 = sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==1));
        yes2fa = sum(p.answer(p.nTrials/2+1:end).*(p.tType(p.nTrials/2+1:end)==2));
        nMisMatching1 = sum(p.tType(1:p.nTrials/2)==1);
        nMatching1 = sum(p.tType(1:p.nTrials/2)==2);
        nMisMatching2 = sum(p.tType(p.nTrials/2+1:end)==1);
        nMatching2 = sum(p.tType(p.nTrials/2+1:end)==2);
        
        hitRate_first(rat,session) = yes1/nMisMatching1;  % a 'yes' (mismatch judgement) on trials that were mismatches (ie, p.tType==1)
        FARate_first(rat,session) = yes1fa/nMatching1; % a 'yes' on matching trials
        
        hitRate_second(rat,session) = yes2/nMisMatching2; % a 'yes' (mismatch judgement) on trials that were mismatches (ie, p.tType==1)
        FARate_second(rat,session) = yes2fa/nMatching2; % a 'yes' on matching trials
                
        
        
        % adjust by adding to only trials
        % THIS IS THE ADJUSTMENT FROM BARENSE ET AL.
        N = 36;
        hitRate_first_adj_some(rat,session) = hitRate_first(rat,session);%#ok<AGROW>
        if hitRate_first_adj_some(rat,session) == 0
            hitRate_first_adj_some(rat,session) = 1/(2*N); %#ok<AGROW>
        elseif hitRate_first_adj_some(rat,session) == 1
            hitRate_first_adj_some(rat,session) = 1 - 1/(2*N); %#ok<AGROW>
        end
        
        FARate_first_adj_some(rat,session) = FARate_first(rat,session); %#ok<AGROW>        
        if FARate_first_adj_some(rat,session) == 0
            FARate_first_adj_some(rat,session) = 1/(2*N); %#ok<AGROW>
        elseif FARate_first_adj_some(rat,session) == 1
            FARate_first_adj_some(rat,session) = 1 - 1/(2*N); %#ok<AGROW>
        end
        
        
        hitRate_second_adj_some(rat,session) = hitRate_second(rat,session); %#ok<AGROW>
        if hitRate_second_adj_some(rat,session) == 0
            hitRate_second_adj_some(rat,session) = 1/(2*N); %#ok<AGROW>
        elseif hitRate_second_adj_some(rat,session) == 1
            hitRate_second_adj_some(rat,session) = 1 - 1/(2*N); %#ok<AGROW>
        end
        
        
        FARate_second_adj_some(rat,session) = FARate_second(rat,session); %#ok<AGROW>  
        if FARate_second_adj_some(rat,session) == 0
            FARate_second_adj_some(rat,session) = 1/(2*N); %#ok<AGROW>
        elseif FARate_second_adj_some(rat,session) == 1
            FARate_second_adj_some(rat,session) = 1 - 1/(2*N); %#ok<AGROW>
        end

    end
end


%%

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


%% Now, actually calculate dPrime as needed
% dPrime_first_raw and dPrime_second_raw => (caudal_LA caudal_HA, PRC_LA, PRC_HA)
close all;

dPrimeTallied(1,1:4) = dPrime_first_adj_some_mean;
dPrimeTallied(2,1:4) = dPrime_second_adj_some_mean;

dPrimePredictions = dPrimeTallied(2,:) - dPrimeTallied(1,:);

barvalues = [dPrime_first_adj_some_mean(4), dPrime_second_adj_some_mean(4) ; dPrime_first_adj_some_mean(3), dPrime_second_adj_some_mean(3)];
errors = [dPrime_first_adj_some_err(4), dPrime_second_adj_some_err(4) ; dPrime_first_adj_some_err(3), dPrime_second_adj_some_err(3)];

figs(23) = figure;
handles = barweb(barvalues...
    , errors...
    , [] ...
    , {'High', 'Low'}...
    , {'Simulation 1'}...
    , {'Stimulus Ambiguity'}...
    , {'d'''}...
    , [rgb('Chocolate') ; rgb('Goldenrod')]);
legend('First Half', 'Second Half','Lesion', 'Location', 'NorthWest');
legend BOXOFF
set(gca,'fontsize',30)
figs(23).CurrentAxes.YLim = [0, 6];
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')

x1 = handles.bars(1).XOffset;
x2 = handles.bars(2).XOffset;

barvalues = [dPrime_first_adj_some_mean(2), dPrime_second_adj_some_mean(2) ; dPrime_first_adj_some_mean(1), dPrime_second_adj_some_mean(1)];
errors = [dPrime_first_adj_some_err(2), dPrime_second_adj_some_err(2) ; dPrime_first_adj_some_err(1), dPrime_second_adj_some_err(1)];

hold on
errorbar([1+x1, 1+x2], barvalues(1,:), errors(1,:), '-kx', 'MarkerSize', 10,'linewidth', 2);
errorbar([2+x1, 2+x2], barvalues(2,:), errors(2,:), '-kx', 'MarkerSize', 10,'linewidth', 2);

saveas(figs(23),[saveFolder, '/dPrime_adj_some'], 'fig');
saveas(figs(23),[saveFolder, '/dPrime_adj_some'], 'jpg');



end

