function [tagName] = mgConcToFdTag(mgConc)
% MGCONCTOFDTAG Returns an FdData tag to label a particular magnesium concentration.
%
% The force-extension data used in the paper is contained in a series of
% FdData objects that all have their magnesium concentration encoded in
% the object's "tags" property. This function converts a magnesium concentration
% to such a tag name.
%
% SYNTAX:
% tagName = mgConcToFdTag(mgConc);
%
% INPUT:
% mgConc = magnesium concentration (in mM; range 0-999).
%
% OUTPUT:
% tagName = FdData tag name.
%
% SEE ALSO:
% fdTagToMgConc, FdData

if (mgConc < 0) || (mgConc > 999)
    error('Invalid argument: "mgConc" out of range.');
end

if mgConc == 0
    tagName = 'buffer';
else
    tagName = sprintf('mg%03d', mgConc);
end

end
