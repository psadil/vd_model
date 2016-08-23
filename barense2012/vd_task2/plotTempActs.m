function [  ] = plotTempActs( prevStimActs, initial_acts, layer,trial,whichCaudal, tType )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
trial
whichCaudal
tType
scrsz = get(groot, 'ScreenSize');

for grid = 1:4
%     figure
    figure('outerposition', scrsz); hold on;
    subplot(1,2,1);
    surf(squeeze(prevStimActs(1,:,:,grid)));
    title('Previous Stim Acts after update')
    subplot(1,2,2);
    surf(squeeze(initial_acts(1,:,:,grid)));
    title('Current Stim Acts before update, Caudal')

end
if layer == 2
    figure; hold on;
    subplot(1,2,1);
    surf(squeeze(prevStimActs(layer,:,:,1)));
    title('Previous Stim Acts after update')
    subplot(1,2,2);
    surf(squeeze(initial_acts(layer,:,:,1)));
    title('Current Stim Acts before update, PRC')
end

close all
end

