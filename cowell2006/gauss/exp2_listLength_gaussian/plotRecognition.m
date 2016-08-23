function [] = plotRecognition(firstRat, lastRat, folderName)
% analyze recognition

% input:
%    firstRat: first rat that was run (usually 1)
%    lastRat: last rat that was run. The range of these two rats are what
%       will be plotted
%   folderName: name of folder that the sessions for these rats can be
%      found in. character string, same as stored in p.nameOfFolder

% output:
%   NA -- but a suite of graphs relating to summarizing the run of these
%      networks that occured in folderName


%%
saveFolder = [pwd,'/graphsAndSession/', folderName];

% load sample data file to get size of parameters
fileName = [saveFolder, '/Session', num2str(1), '_Rat', num2str(1)];
load(fileName)

numRats = lastRat-firstRat+1;

recog = zeros(numRats,p.nSess/2,2,max(p.nTrials),p.nStimSets);
recog_gauss = zeros(size(recog));
recogByLayer = zeros(numRats,p.nSess,2,max(p.nTrials),p.nStimSets);
recogByLayer_gauss = zeros(size(recogByLayer));

for rat = firstRat:lastRat
    for session = 1:p.nSess
        
        fileName = [saveFolder, '/Session', num2str(session), '_Rat', num2str(rat)];
        load(fileName)
        fprintf ('\nloading rat %d, session %d', rat, session);
        
        
        %------------------------------------------------------------------
        % collect all trial-wise selectivity of choice phase
        %------------------------------------------------------------------
        
        recogByLayer(rat,session,1,1:p.nTrials(p.stimCond),:) = p.recogByLayer(:,1,:);
        recogByLayer_gauss(rat,session,1,1:p.nTrials(p.stimCond),:) = p.recogByLayer_gauss(:,1,:);
        if p.layer == 1
            recog(rat,p.stimCond,1,1:p.nTrials(p.stimCond),:) = p.recognition;
            recog_gauss(rat,p.stimCond,1,1:p.nTrials(p.stimCond),:) = p.recognition_gauss;
        else
            recogByLayer(rat,session,2,1:p.nTrials(p.stimCond),:) = p.recogByLayer(:,2,:);
            recogByLayer_gauss(rat,session,2,1:p.nTrials(p.stimCond),:) = p.recogByLayer_gauss(:,2,:);
            recog(rat,p.stimCond,2,1:p.nTrials(p.stimCond),:) = p.recognition ;
            recog_gauss(rat,p.stimCond,2,1:p.nTrials(p.stimCond),:) = p.recognition_gauss;
        end
    end
end

%%

% first, find the average for a rat in a given condition
rats = zeros(numRats,length(p.nTrials),2);
for stimCond = 1:length(p.nTrials)
    rats(:,stimCond,:) = squeeze(mean(squeeze(mean(recog(:,stimCond,:,1:p.nTrials(stimCond),:),4)),3));
    
end
recog_mean = squeeze(mean(rats,1));
recog_sem = squeeze(std(rats,1) ./ sqrt(numRats));


% gauss recognition
rats_gauss = zeros(numRats,length(p.nTrials),2);
for stimCond = 1:length(p.nTrials)
    rats_gauss(:,stimCond,:) = squeeze(mean(squeeze(mean(recog_gauss(:,stimCond,:,1:p.nTrials(stimCond),:),4)),3));
end
recog_gauss_mean = squeeze(mean(rats_gauss,1));

recog_gauss_sem = squeeze(std(rats,1) ./ sqrt(numRats));

% -------------------------------------------------------------------------
% recognition by layer
% -------------------------------------------------------------------------

% sigmoid
rats_byLayer = zeros(numRats,p.nSess,2);
for sess = 1:p.nSess
    
    % fill up with recognition calculated from each layer. Will be some 0s,
    % since PRC was only available on sessions 5-8
    if sess <= length(p.nTrials)
        rats_byLayer(:,sess,1) = squeeze(mean(recogByLayer(:,sess,1,1:p.nTrials(sess)),4));
    elseif sess > length(p.nTrials)
        idx = sess - length(p.nTrials);
        rats_byLayer(:,sess,:) = squeeze(mean(squeeze(mean(recogByLayer(:,sess,:,1:p.nTrials(idx),:),4)),3));
    end
end
recogByLayer_mean = squeeze(mean(rats_byLayer,1));
recogByLayer_sem = squeeze(std(rats_byLayer,1) ./ sqrt(numRats));


% gauss
rats_byLayer_gauss = zeros(numRats,p.nSess,2);
for sess = 1:p.nSess
    
    % fill up with recognition calculated from each layer. Will be some 0s,
    % since PRC was only available on sessions 5-8
    if sess <= length(p.nTrials)
        rats_byLayer_gauss(:,sess,1) = squeeze(mean(recogByLayer_gauss(:,sess,1,1:p.nTrials(sess)),4));
    elseif sess > length(p.nTrials)
        idx = sess - length(p.nTrials);
        rats_byLayer_gauss(:,sess,:) = squeeze(mean(squeeze(mean(recogByLayer_gauss(:,sess,:,1:p.nTrials(idx),:),4)),3));
    end
end
recogByLayer_gauss_mean = squeeze(mean(rats_byLayer_gauss,1));
recogByLayer_gauss_sem = squeeze(std(rats_byLayer_gauss,1) ./ sqrt(numRats));


%%
close all

figs(1) = figure;
hold on
plot(1:4,recog_mean(:,1), '--ok', 'MarkerSize',10)
plot(1:4,recog_mean(:,2),'-ok','MarkerFaceColor','k', 'MarkerSize',10)
ax = gca;
ax.XTick = 1:4;
legend('Lesion','Control')
legend('boxoff')

saveas(figs(1),[saveFolder, '/recog'],'fig');
saveas(figs(1),[saveFolder, '/recog'],'jpg');

% gauss activation
figs(2) = figure;
hold on
plot(1:4,recog_gauss_mean(:,1), '--ok', 'MarkerSize',10)
plot(1:4,recog_gauss_mean(:,2),'-ok','MarkerFaceColor','k', 'MarkerSize',10)
ax = gca;
ax.XTick = 1:4;
legend('Lesion','Control')
legend('boxoff')

saveas(figs(2),[saveFolder, '/recogGauss'],'fig');
saveas(figs(2),[saveFolder, '/recogGauss'],'jpg');



% recognition, sigmoidal
figs(3) = figure;
subplot(1,2,1)
barweb(recog_mean(:,2), recog_sem(:,2), [], {'control'})
xlabel('stim condition');
ylabel('recognition sigmoid');
legend({'1','6','12','18'},'Location','best');
figs(3).CurrentAxes.YLim = [0,...
    max(recog_mean(:))+max(recog_sem(:))];

subplot(1,2,2)
barweb(recog_mean(:,1), recog_sem(:,1), [], {'lesion'})
xlabel('stim condition');
ylabel('recognition sigmoid');
legend({'1','6','12','18'},'Location','best');
figs(3).CurrentAxes.YLim = [0,...
    max(recog_mean(:))+max(recog_sem(:))];

saveas(figs(3),[saveFolder, '/recog_sigm_barweb'],'fig');
saveas(figs(3),[saveFolder, '/recog_sigm_barweb'],'jpg');


% recognition, gaussian
figs(4) = figure;
subplot(1,2,1)
barweb(recog_gauss_mean(:,2), recog_gauss_sem(:,2), [], {'control'})
xlabel('stim condition');
ylabel('recognition gauss');
legend({'1','6','12','18'},'Location','best');
figs(4).CurrentAxes.YLim = [0,...
    max(recog_gauss_mean(:))+max(recog_gauss_sem(:))];

subplot(1,2,2)
barweb(recog_gauss_mean(:,1), recog_gauss_sem(:,1), [], {'lesion'})
xlabel('stim condition');
ylabel('recognition gauss');
legend({'1','6','12','18'},'Location','best');
figs(4).CurrentAxes.YLim = [0,...
    max(recog_gauss_mean(:))+max(recog_gauss_sem(:))];


saveas(figs(4),[saveFolder, '/recog_gauss_barweb'],'fig');
saveas(figs(4),[saveFolder, '/recog_gauss_barweb'],'jpg');


%% recognition by layer

% sigmoid
figs(5) = figure;
subplot(1,2,1)
barweb([recogByLayer_mean(1:4,2),recogByLayer_mean(5:8,2)]',...
    [recogByLayer_sem(1:4,2),recogByLayer_sem(5:8,2)]', [], {'PRC, lesion', 'PRC, control'})
xlabel('stim condition');
ylabel('recognition sigmoid');
legend({'1','6','12','18'},'Location','best');
figs(5).CurrentAxes.YLim = [0,...
    max(recogByLayer_mean(:,2))+max(recogByLayer_sem(:,2))];

subplot(1,2,2)
barweb([recogByLayer_mean(1:4,1),recogByLayer_mean(5:8,1)]',...
    [recogByLayer_sem(1:4,1),recogByLayer_sem(5:8,1)]', [], {'caudal, lesion', 'caudal, control'})
xlabel('stim condition');
ylabel('recognition sigmoid');
legend({'1','6','12','18'},'Location','best');
figs(5).CurrentAxes.YLim = [0,...
    max(recogByLayer_mean(:,1))+max(recogByLayer_sem(:,1))];

saveas(figs(5),[saveFolder, '/recogByLayer_sigm_barweb'],'fig');
saveas(figs(5),[saveFolder, '/recogByLayer_sigm_barweb'],'jpg');

% gaussian
figs(6) = figure;
subplot(1,2,1)
barweb([recogByLayer_gauss_mean(1:4,2),recogByLayer_gauss_mean(5:8,2)]',...
    [recogByLayer_gauss_sem(1:4,2),recogByLayer_gauss_sem(5:8,2)]',...
    [], {'PRC, lesion', 'PRC, control'})
xlabel('stim condition');
ylabel('recognition gaussina');
legend({'1','6','12','18'},'Location','best');
figs(6).CurrentAxes.YLim = [0,...
    max(recogByLayer_gauss_mean(:,2))+max(recogByLayer_gauss_sem(:,2))];

subplot(1,2,2)
barweb([recogByLayer_gauss_mean(1:4,1),recogByLayer_gauss_mean(5:8,1)]',...
    [recogByLayer_gauss_sem(1:4,1),recogByLayer_gauss_sem(5:8,1)]',...
    [], {'caudal, lesion', 'caudal, control'})
xlabel('stim condition');
ylabel('recognition gaussian');
legend({'1','6','12','18'},'Location','best');
figs(6).CurrentAxes.YLim = [0,...
    max(recogByLayer_gauss_mean(:,1))+max(recogByLayer_gauss_sem(:,1))];

saveas(figs(6),[saveFolder, '/recogByLayer_gauss_barweb'],'fig');
saveas(figs(6),[saveFolder, '/recogByLayer_gauss_barweb'],'jpg');


end

