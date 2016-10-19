function [F] = ftWLCInv_f0(d, Lp, Lc, S, C, g0, g1, Fc, F0, kT)
% ftWLCInv, with force offset included (for convenience)
%
% SEE ALSO:
% ftWLCInv

if nargin < 10
    kT = 4.11;
end

F = ftWLCInv(d, Lp, Lc, S, C, g0, g1, Fc, kT)+F0;

end


