function plotfdfits(fdc, fits, model)
% PLOTFDFITS Plot fits to a FdDataCollection
%
% SYNTAX:
% plotfdfits(fdc, fits, model)
%
% INPUT:
% fdc = the FdDataCollection that was passed into 'fitfd'
% fits = the output from 'fitfd'
% model = the model object used for fitting

explorecell(fdc.items, @(ax,fd,idx) plotfdfit(ax, fd, fits(idx).fitobject, model));

end
