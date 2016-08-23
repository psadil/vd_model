function [inp_matrix] = VD_gen_limited_input(numInputDims,p)
numRows = p.numRows;

% avail_features = [ .1 .1 .1; .1 .1 .9; ...
%                    .1 .9 .1; .1 .9 .9; ...
%                    .9 .1 .1; .9 .1 .9; ...
%                    .9 .9 .1; .9 .9 .9];
%
nSimpleConj = (p.numInputDims_Caudal / p.nDimReps) ^ p.nStimFactors;
count=1;
avail_features = zeros(nSimpleConj, p.numInputDims_Caudal / p.nDimReps);
for inp1 = 1:p.numGrids_Caudal,
    for inp2 = 1:p.nStimFactors,
        avail_features(count,:) = [inp1 inp2];
        count=count+1;
    end
end

% count=1;
% avail_features = zeros(p.components, p.numInputDims_Caudal);
% % for inp1 = 1:nSimpleFeats,
% for inp2 = 1:p.numInputDims_Caudal,
%     %         availFeat(count,:) = [inp1 inp2];
%     avail_features(count) = inp2;
%     count=count+1;
%     %     end
% end


avail_features(avail_features==1)=0.05;
avail_features(avail_features==2)=0.35;
avail_features(avail_features==3)=0.65;
avail_features(avail_features==4)=0.95;


%%%%% Generate an input
% if numInputDims == 1,
%     inp_matrix = avail_features(ceil(rand*8),1)*ones(numRows,NUM_COLS);
%   %  inp_vector = inp_matrix(1,1,1);
% elseif numInputDims == 2,
%     segment = avail_features(ceil(rand*8),:);
%     inp_matrix = cat(3, segment(1)*ones(numRows,NUM_COLS), segment(2)*ones(numRows));
% 	      % Take the input x,y coords, and make each fill a slice
%               % of a matrix the same size as m, so that they can be subtracted
%  %   inp_vector = [inp_matrix(1,1,1) inp_matrix(1,1,2)];
if numInputDims == p.numInputDims_Caudal / p.nDimReps,
    segment = avail_features(randi(nSimpleConj),:);
    
    %     segment = zeros(p.numInputDims_Caudal);
    %     for seg = 1:p.numInputDims_Caudal
    %         segment(seg,:) = avail_features(randi(p.components),:);
    %     end
    
    inp_matrix_temp = ones(numRows,numRows,p.numInputDims_Caudal / p.nDimReps);
    inp_matrix = zeros(size(inp_matrix_temp));
    for dim = 1:p.numInputDims_Caudal/p.nDimReps
        inp_matrix(:,:,dim) = segment(dim) * inp_matrix_temp(:,:,dim);
    end
    
    
    
    %     inp_matrix = cat(3, segment(1)*ones(numRows), segment(2)*ones(numRows), segment(3)*ones(numRows));
    % Take the input x,y coords, and make each fill a slice
    % of a matrix the same size as m, so that they can be subtracted
    %  inp_vector = [inp_matrix(1,1,1) inp_matrix(1,1,2) inp_matrix(1,1,3) inp_matrix(1,1,4) inp_matrix(1,1,5) inp_matrix(1,1,6)];
elseif numInputDims == p.numInputDims_PRC / p.nDimReps
    for seg = 1:p.numGrids_Caudal ,
        segment(seg,:) = avail_features(randi(nSimpleConj),:);
    end
    
    inp_matrix_temp = ones(numRows,numRows,p.numInputDims_PRC / p.nDimReps);
    inp_matrix = zeros(size(inp_matrix_temp));
    for dim = 1:p.numInputDims_PRC / p.nDimReps
        inp_matrix(:,:,dim) = segment(dim) * inp_matrix_temp(:,:,dim);
    end
    
    %     inp_matrix = cat(3, segment(1,1)*ones(numRows), segment(1,2)*ones(numRows), segment(1,3)*ones(numRows), ...
    %                         segment(2,1)*ones(numRows), segment(2,2)*ones(numRows), segment(2,3)*ones(numRows), ...
    %                         segment(3,1)*ones(numRows), segment(3,2)*ones(numRows), segment(3,3)*ones(numRows), ...
    %                         segment(4,1)*ones(numRows), segment(4,2)*ones(numRows), segment(4,3)*ones(numRows), ...
    %                         segment(5,1)*ones(numRows), segment(5,2)*ones(numRows), segment(5,3)*ones(numRows));
else
    fprintf('Invalid number of stimulus dimensions \n');
end

inp_matrix = reshape(repmat(reshape(inp_matrix,...
    [1,p.numRows,1,p.numRows,1,numInputDims]),[1,1,1,1,p.nDimReps]),...
    [p.numRows,p.numRows,numInputDims*p.nDimReps]);


end