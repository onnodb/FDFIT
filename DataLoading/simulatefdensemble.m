function [fdc] = simulatefdensemble(N, varargin)
% SIMULATEFDENSEMBLE Generate an ensemble of simulated F,d data
%
% This is a simple wrapper around 'simulatefddata', which is
% simply called N times. All arguments are passed directly
% in to 'simulatefddata'.
%
% NOTE: Uses a 'parfor' loop to parallelize the process and increase performance.
%
% SYNTAX
% fdc = simulatefdensemble(..., 'key', 'value', ...)
%
% SEE ALSO:
% simulatefddata

items = cell(N,1);
parfor i = 1:N
    items{i} = simulatefddata(varargin{:}); %#ok
end

fdc = FdDataCollection('skipduplicatecheck', items);

end
