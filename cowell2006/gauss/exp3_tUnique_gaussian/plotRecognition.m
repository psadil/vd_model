function [] = plotRecognition(firstRat, lastRat, folderName)
% analyze recognition

% scrsz = get(groot, 'ScreenSize');

saveFolder = [pwd,'/graphsAndSession/', folderName];

% load sample data file to get size of parameters
fileName = [saveFolder, '/Session', num2str(1), '_Rat', num2str(1)];
load(fileName)

numRats = lastRat-firstRat+1;

recognition = zeros(numRats,p.nSess/2,2,p.nTrials);

for rat = firstRat:lastRat
    for session = 1:p.nSess
        
        fileName = [saveFolder, '/Session', num2str(session), '_Rat', num2str(rat)];
        load(fileName)
        fprintf ('\nloading rat %d, session %d', rat, session);
        
        
        %------------------------------------------------------------------
        % collect all trial-wise selectivity of choice phase
        %------------------------------------------------------------------
        
        if session <= p.nSess/2 % lesion trials
            recognition(rat,session,1,:) = p.recognition;
        else % control trials
            recognition(rat,session-p.nSess/2,2,:) = p.recognition;
        end
    end
end

%%

% first average across trial, then across rats
meanRecog = squeeze(mean(squeeze(mean(recognition,4)),1));

% take mean of 4 trials (4th dim), then std across rats, and divide by
% sqrt(numRats)
recog_sem = squeeze(std(squeeze(mean(recognition,4)),1,1))./sqrt(numRats);


% convert to table for writting to csv
% recogT = table( squeeze(mean(recognition),1));
% recogT.Properties.VariableNames

% write that table as csv

%%
close all



figs(4) = figure;
barweb([meanRecog(:,2),meanRecog(:,1)], [recog_sem(:,2),recog_sem(:,1)], [], {'trial unique','repeating'})
xlabel('stim condition');
ylabel('recognition');
legend({'control', 'lesion'},'Location','best');
figs(4).CurrentAxes.YLim = [min(meanRecog(:,2))-max(recog_sem(:,2)),...
    max(meanRecog(:,2))+max(recog_sem(:,2))];

saveas(figs(4),[saveFolder, '/recog_barweb'],'fig');
saveas(figs(4),[saveFolder, '/recog_barweb'],'jpg');


end

