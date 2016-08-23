function [ p ] = createStimOrder( p )
%VD_createStimOder creates stimulus order to be presented to simulation
%during a given session
%   called during run_sim



p.tType = cat(2,ones(1,p.nMismatch), 2*ones(1,p.nMatch));
p.tType = p.tType(randperm(length(p.tType)));
%check to see whether there are more than 3 trials in a row the same % if yes, redo
check = 1;
while check
    onesInSeries = p.tType==1;
    fourInaRow = onesInSeries(1:end-3)+onesInSeries(2:end-2)+onesInSeries(3:end-1)+onesInSeries(4:end);
    anyFours=fourInaRow==4;
    sum(anyFours);
    twosInSeries = p.tType==2;
    fourTwosInaRow = twosInSeries(1:end-3)+twosInSeries(2:end-2)+twosInSeries(3:end-1)+twosInSeries(4:end);
    anyFours2=fourTwosInaRow==4;
    sum(anyFours2);
    if sum(anyFours)>0 || sum(anyFours2)>0
        p.tType = p.tType(randperm(length(p.tType)));
        check=1;
    else
        check=0;
    end
end
p.stimOrder = [randperm(p.nMismatch)' randperm(p.nMatch)'];




end

