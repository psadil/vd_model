function [inp_mat] = gen_limited_input(numInputDims,numRows)
% generate and normalize random input or weights


inp = rand(numRows,numRows,numInputDims);

inp_mat = inp ./...
    repmat(squeeze(sqrt(sum(inp.^2,3))),[1,1,numInputDims]);

% to test, dot product of each node with itself (eg. n1 = squeeze(inp(1,8,:))) should
% be 1; dot(n1,n1) == 1;

end