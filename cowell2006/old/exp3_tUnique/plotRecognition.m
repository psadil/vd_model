function [] = plotRecognition(firstRat, lastRat, folderName)
% analyze recognition

% scrsz = get(groot, 'ScreenSize');


%%
saveFolder = [pwd,'/graphsAndSession/', folderName];

% load sample data file to get size of parameters
fileName = [saveFolder, '/Session', num2str(1), '_Rat', num2str(1)];
load(fileName)

numRats = lastRat-firstRat+1;

recog = zeros(numRats,p.nSess/2,2,p.nTrials/2);
recogByLayer = zeros(numRats,p.nSess,2,p.nTrials/2);


%%
for rat = firstRat:lastRat
    for session = 1:p.nSess
        
        fileName = [saveFolder, '/Session', num2str(session), '_Rat', num2str(rat)];
        load(fileName)
        fprintf ('\nloading rat %d, session %d', rat, session);
        
        
        %------------------------------------------------------------------
        % collect all trial-wise selectivity of choice phase
        %------------------------------------------------------------------
        
        recogByLayer(rat,session,1,:) = p.recogByLayer(:,1);
        if p.layer == 1
            recog(rat,p.stimCond,1,:) = p.recognition;
        else
            recogByLayer(rat,session,2,:) = p.recogByLayer(:,2);
            recog(rat,p.stimCond,2,:) = p.recognition;
        end
    end
end

%%

% first average across trial, then across rats
% gives stimCond x layer
recog_mean = squeeze(mean(squeeze(mean(recog,4)),1));

% take mean of 4 trials (4th dim), then std across rats, and divide by
% sqrt(numRats)
recog_sem = squeeze(std(squeeze(mean(recog,4)),1))./sqrt(numRats);

% -------------------------------------------------------------------------
% by layer
% -------------------------------------------------------------------------

recogByLayer_mean = squeeze(mean(squeeze(mean(recogByLayer,4)),1));

% take mean of 4 trials (4th dim), then std across rats, and divide by
% sqrt(numRats)
recogByLayer_sem = squeeze(std(squeeze(mean(recogByLayer,4)),1))./sqrt(numRats);


%%
close all

figs(1) = figure;
hold on
plot(1:2,recog_mean(:,1), '--ok', 'MarkerSize',10)
plot(1:2,recog_mean(:,2),'-ok','MarkerFaceColor','k', 'MarkerSize',10)
ax = gca;
ax.XTick = 1:2;
legend('Lesion','Control')
legend('boxoff')

saveas(figs(1),[saveFolder, '/recog'],'fig');
saveas(figs(1),[saveFolder, '/recog'],'jpg');

% barplot, to put error bars easily
figs(2) = figure;
barweb([recog_mean(:,2),recog_mean(:,1)],...
    [recog_sem(:,2),recog_sem(:,1)], [], {'trial unique','repeating'})
xlabel('stim condition');
ylabel('recognition');
legend({'control', 'lesion'},'Location','best');
figs(2).CurrentAxes.YLim = [min(recog_mean(:))-max(recog_sem(:)),...
    max(recog_mean(:))+max(recog_sem(:))];

saveas(figs(2),[saveFolder, '/recog_barweb'],'fig');
saveas(figs(2),[saveFolder, '/recog_barweb'],'jpg');

% by layer
figs(3) = figure;
subplot(1,2,1)
barweb([recogByLayer_mean(1:2,2),recogByLayer_mean(3:4,2)]',...
    [recogByLayer_sem(1:2,2),recogByLayer_sem(3:4,2)]', [], {'PRC, lesion', 'PRC, control'})
xlabel('stim condition');
ylabel('recognition sigmoid');
legend({'trial unique', 'repitition'},'Location','best');
figs(3).CurrentAxes.YLim = [min(recog_mean(:))-max(recog_sem(:)),...
    max(recogByLayer_mean(:))+max(recogByLayer_sem(:))];

subplot(1,2,2)
barweb([recogByLayer_mean(1:2,1),recogByLayer_mean(3:4,1)]',...
    [recogByLayer_sem(1:2,1),recogByLayer_sem(3:4,1)]', [], {'caudal, lesion', 'caudal, control'})
xlabel('stim condition');
ylabel('recognition sigmoid');
legend({'trial unique', 'repitition'},'Location','best');
figs(3).CurrentAxes.YLim = [min(recog_mean(:))-max(recog_sem(:)),...
    max(recogByLayer_mean(:))+max(recogByLayer_sem(:))];

saveas(figs(3),[saveFolder, '/recogByLayer_barweb'],'fig');
saveas(figs(3),[saveFolder, '/recogByLayer_barweb'],'jpg');

end