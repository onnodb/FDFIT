function [d] = fOdijk_d0(F, Lp, Lc, S, d0, kT)

if nargin < 6
    kT = 4.11;
end

d = fOdijk(F, Lp, Lc, S, kT)-d0;

end
