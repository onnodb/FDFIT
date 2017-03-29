function plotglobalfdfits(fdc, fitRes, model)
% PLOTGLOBALFDFITS Plot globally fitted force-extension data with associated fits.
%
% SYNTAX:
% plotglobalfdfits(fdc, fitRes)
%
% INPUT:
% fdc = the FdDataCollection the global fit was performed on.
% fitRes = the output of the "fitfdglobal" function.

explorecell(fdc.items, @plotit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [] = plotit(ax, fd, idx)
        fn = fieldnames(fitRes);
        p = struct();
        for i = 1:length(fn)
            if isscalar(fitRes.(fn{i}))
                p.(fn{i}) = fitRes.(fn{i});
            else
                p.(fn{i}) = fitRes.(fn{i})(idx);
            end
        end

        plotfdfit(ax, fd, p, model);
    end


end
