function [fdc_aligned] = alignFdCurves(fdc, varargin)
% ALIGNFDCURVES Aligns a set of F,d curves, prior to a TWLC analysis.
%
% Uses a global fit of the eWLC ("odijk-d0-f0" model) to the given set of
% F,d curves to align them, both along the "d" and the "F" axes.
%
% First, the global fit is used to correct distance and force offsets. Then,
% (optionally), curves are scaled so that their overstretching regimes
% overlap (using the average force value in this region for all curves
% as a target).
%
% SYNTAX:
% fdc_aligned = alignFdCurves(fdc);
%
% INPUT:
% fdc = an FdDataCollection with F,d curves to align
%
% KEY-VALUE PAIR ARGUMENTS:
% fitfdglobalArgs = optional cell array with arguments to pass into the
%       'fitfdglobal' function (such as a value for 'Lc').
% OSPlateauRegion = data region along the distance axis [min max], used for
%       scaling the curves vertically (along the F axis) so that their
%       overstretching regimes overlap. If not given, this step is not
%       performed. (Default: empty = skip this step)
%
% OUTPUT:
% fdc_aligned = FdDataCollection with aligned F,d curves
%
% SEE ALSO:
% fitfdglobal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse & validate input

defArgs = struct(...
                  'fitfdglobalArgs',                        {{}} ...
                , 'OSPlateauRegion',                        [] ...
                );
args = parseArgs(varargin, defArgs);

if ~isempty(args.OSPlateauRegion) && (~isnumeric(args.OSPlateauRegion) || ...
        ~isvector(args.OSPlateauRegion) || (length(args.OSPlateauRegion) ~= 2))
    error('Invalid arguments "OSPlateauRegion": length-2 vector [min max] expected');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shift along along "d" & "F" axes for baseline alignment

% Perform a global fit of the "odijk-d0-f0" model ("Odijk" model with fixed
% Lc, with an additional distance offset and force offset).
globalFit = fitfdglobal(fdc, args.fitfdglobalArgs{:});

% Loop over global fit results, and correct found offsets for all F,d curves.
fdc_shifted = FdDataCollection();
for i = 1:fdc.length
    fdc_shifted.add('skipduplicatecheck', ...
                    fdc.items{i}.shift('d', globalFit.d0(i)).shift('f', -globalFit.F0(i)) ...
                    );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scale to match overstretching plateaus (if requested)

if isempty(args.OSPlateauRegion)
    fdc_aligned = fdc_shifted;
else
    % Calculate mean plateau height within given plateau region for each curve.
    plateauHeight = zeros(fdc.length,1);
    for i = 1:fdc_shifted.length
        cur_f = fdc_shifted.items{i}.f;
        cur_d = fdc_shifted.items{i}.d;

        plateauHeight(i) = mean(cur_f(cur_d >= args.OSPlateauRegion(1) & cur_d <= args.OSPlateauRegion(2)));

        if isnan(plateauHeight(i))
            % NaN value for mean: this means there was no data for the OS
            % region in this dataset. That's an error!
            error('No OS plateau data for "%s" (#%d)', fdc.items{i}.name, i);
        end
    end

    meanPlateauHeight = mean(plateauHeight);

    % Now rescale each curve to match the calculated plateau height
    fprintf('Stretching F,d curves to match average plateau height: %g\n', meanPlateauHeight);

    fdc_aligned = FdDataCollection();
    for i = 1:fdc_shifted.length
        fdc_aligned.add('skipduplicatecheck', ...
                        fdc_shifted.items{i}.scale('f', meanPlateauHeight/plateauHeight(i)) ...
                        );
    end

    % Oops, this rescaling step messed up the baseline alignment :)
    % So do another round of that (this time *without* scaling) to fix it.
    fdc_aligned = alignFdCurves(fdc_aligned, 'fitfdglobalArgs', args.fitfdglobalArgs);
end

end

