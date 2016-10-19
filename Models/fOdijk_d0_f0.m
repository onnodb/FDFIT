function [d] = fOdijk_d0_f0(F, Lp, Lc, S, d0, F0, kT)

if nargin < 7
    kT = 4.11;
end

d = fOdijk(F-F0, Lp, Lc, S, kT)-d0;

end
