classdef TwlcAnalysis < handle
    % TWLCANALYSIS Encapsulation of the analysis of the tWLC.
    %
    % This is a convenience class that encapsulates the analysis procedure
    % for the paper: extracting tWLC fit parameters for a dataset.

    properties

        % ----- Input data

        % FdDataCollection containing the data to be analyzed.
        % Each FdData object should contain a clean F,d curve.
        data                                    = [];


        % ----- Cache

        % Aligned versions of the input data. Will be generated
        % automatically if not given. If you have already performed
        % alignment yourself, you can also pass in this aligned version
        % by using this property.
        dataAligned                             = [];


        % ----- Sweeping & bootstrapping settings

        % Contour length of the DNA (um).
        Lc                                      = 16.5;

        % tWLC critical force (default: 30.6 pN).
        % This is the force at which the twist-stretch coupling changes from
        % a constant to a linear function (g0 + g1 F).
        Fc                                      = 30.6;

        % Twist rigidity (default: 440 pN nm^2).
        C                                       = 440;

        % Number of bootstrap iterations to perform for error bars.
        nBootstrapIter                          = 100;

        % Number of sweeps to perform to determine Fmax.
        nSweeps                                 = 100;

        % Vector [Fmax_min Fmax_max]. "nSweeps" sweeps are performed
        % for Fmax values between "Fmax_min" and "Fmax_max" (pN).
        sweepBoundaries                         = [40 70];

        % Distance region for the overstretching plateau.
        % Should be a vector [Dmin Dmax], with Dmin and Dmax both
        % in um.
        osPlateauRegion                         = [18.5 20];

        % 2x4 matrix [Lp_min S_min g0_min g1_min; Lp_max ...]
        % indicating fit parameter bounds.
        twlcFitBounds                           = [5   500  -2000 2  ; ... % Lp S g0 g1
                                                   150 5000 -100  60];

        % Optional vector with starting values for the fit parameters
        % [Lp S g0 g1].
        twlcFitStartParams                      = [];


        % ----- Misc.

        % Parameter names (for plotting; usually no need to change this).
        params = {  'Lp',       'L_p (nm)',     'nm'        , [] ; ...
                    'S',        'S (pN)',       'pN'        , [] ; ...
                    'g0',       'g0 (pN nm)',   'pN nm'     , [] ; ...
                    'g1',       'g1 (nm)',      'nm'        , []   ...
                    };

    end

    properties (SetAccess = private)

        % Buffer for the sweep results (internal use).
        sweepRes                                = [];

        % Buffer for the sweep results (internal use).
        sweepData                               = [];

    end

    % ------------------------------------------------------------------------

    methods

        function [self] = TwlcAnalysis(varargin)
            % TWLCANALYSIS Constructor.
            %
            % SYNTAX:
            % twa = TwlcAnalysis('key', value, ...)
            %
            % KEY-VALUE PAIR ARGUMENTS:
            % All class properties can be initialized using the key-value pair
            % system. (Note: these are case-sensitive).

            parseClassArgs(varargin, self);

            if isempty(self.data) || ~isa(self.data, 'FdDataCollection') || (self.data.length == 0)
                error('No data or invalid data given: FdDataCollection expected.');
            end
        end

        function analyze(self)
            % ANALYZE Runs the sweeps & bootstraps.
            %
            % Note that you don't need to call this method again if you're just
            % tweaking the cut-off point determination.

            %% Align the F,d curves (all together).
            if isempty(self.dataAligned)
                disp('Aligning...');
                self.dataAligned = alignFdCurves(self.data, 'OSPlateauRegion', self.osPlateauRegion);
            else
                disp('Using pre-set aligned data.');
            end

            [self.sweepRes, self.sweepData] = self.analyzeData(self.dataAligned);

            self.disp();
        end

        function disp(self)
            % DISP Display results.

            disp('<a href="matlab:doc TwlcAnalysis;">TwlcAnalysis</a> object');
            fprintf('tWLC analysis on %d force-extension curves.\n', self.data.length);

            disp('Parameters:');
            for p = {'Lc','Fc','C','nBootstrapIter','nSweeps','sweepBoundaries',...
                     'osPlateauRegion','twlcFitBounds','twlcFitStartParams'}
                fprintf('    %s = %s\n', p{1}, num2str(self.(p{1})));
            end

            if isempty(self.sweepRes)
                disp(' --> Call "analyze" first to run the analysis.');
            else
                disp('Fit results:');
                nParams  = size(self.params,1);

                [cutoffForce, cutoffIdx]  = self.findCutOffForce();
                fprintf('  Fmax: %g pN\n', cutoffForce);
                for iParam = 1:nParams
                    paramVals = squeeze(self.sweepData(iParam,:,cutoffIdx));
                    fprintf('    %2s: %g +/- %g %s\n', ...
                            self.params{iParam,1}, ...
                            mean(paramVals), std(paramVals), ...
                            self.params{iParam,3} ...
                            );
                end

                disp('Call "plotFit" to plot the fit, and "plotSweeps" to plot the');
                disp('sweep results.');
            end
        end

        function [cutoffForce, cutoffIdx] = findCutOffForce(self)
            % FINDCUTOFFFORCE Finds the cut-off force
            %
            % Determines the upper limit of validity of the tWLC.
            %
            % OUTPUT:
            % cutoffForce = upper limit of validity of the tWLC, i.e., the
            %       maximum force that can be used for fitting the tWLC (pN).
            % cutoffIdx = index of the sweep that corresponds to 'cutoffForce'
            %       (index into 'sweepRes.options.sweepPosRight).

            if isempty(self.sweepData)
                error('No analysis data; run "analyze" first.');
            end

            % Cut-off force determined from maximum GOF (goodness-of-fit).
            gofVals = zeros(self.nBootstrapIter, self.nSweeps);
            for i = 1:self.nBootstrapIter
                for j = 1:self.nSweeps
                    gofVals(i,j) = self.sweepRes{i}.params.err(j);
                end
            end
            gofVals = mean(gofVals);
            [~, maxGofIdx] = max(gofVals);
            cutoffForce = self.sweepRes{1}.options.sweepPosRight(maxGofIdx);
            cutoffIdx   = maxGofIdx;
        end

        function plotFit(self)
            % PLOTFIT Plot the fit results.

            if isempty(self.sweepData)
                error('No analysis data; run "analyze" first.');
            end

            allFd = self.dataAligned.concatenatedData();

            [cutoffForce, cutoffIdx] = self.findCutOffForce();
            paramVals = zeros(4,1);
            for iParam = 1:4
                paramVals(iParam) = mean(squeeze(self.sweepData(iParam,:,cutoffIdx)));
            end

            figure;
            plotfd(allFd);
            hold('on');
            dVals = sort(allFd.d);
            plot(...
                dVals, ...
                ftWLCInv(dVals, paramVals(1), self.Lc, paramVals(2), ...
                         self.C, paramVals(3), paramVals(4), self.Fc), ...
                '-r' ...
                );
            plot(...
                dVals, ...
                cutoffForce .* ones(size(dVals)), ...
                '--', 'Color', [0.7 0.7 0.7] ...
                );
        end

        function plotSweeps(self, varargin)
            % PLOTSWEEPS Plot the sweep results.
            %
            % Plots the fit results as a function of the maximum/cut-off force
            % Fmax. The red marker indicates the final cut-off force found,
            % i.e., the value of Fmax that produces the optimal goodness-of-fit
            % (and the value of Fmax shown in the fit results).
            %
            % KEY-VALUE PAIR ARGUMENTS:
            % hFig = optional figure handle; the figure is plotted in this
            %       window.
            % ylims = optional 4x2 matrix, with each row giving a vector
            %       [ymin ymax] to pass into "ylim". The ordering of the rows
            %       is: Lp, S, g0, g1.
            %
            % FLAG ARGUMENTS
            % plotErr = if given, also plots an 'error' graph (using the 'err'
            %       parameter of the sweep results, which usually corresponds
            %       to 'goodness-of-fit'; see 'sweepfdfits', the 'gofType'
            %       argument).

            defArgs = struct(...
                              'hFig',                   [] ...
                            , 'plotErr',                false ...
                            , 'ylims',                  [] ...
                            );
            args = parseArgs(varargin, defArgs, {'plotErr'});

            if isempty(self.sweepRes)
                error('No data to plot; run "analyze" first.');
            end

            if isempty(args.hFig)
                figure;
            else
                figure(args.hFig);
                clf();
            end

            %% Create plot
            nParams = size(self.params,1);
            if args.plotErr
                nParams = nParams + 1;
            end
            nRows = floor(sqrt(nParams));
            nCols = ceil(nParams/nRows);

            for iParam = 1:nParams
                subplot(nRows,nCols,iParam);

                if args.plotErr && (iParam == nParams)
                    makeErrPlot();
                else
                    makeParamPlot(iParam);
                end
            end

            % >> nested function
            function makeParamPlot(paramIdx)
                paramName  = self.params{paramIdx,1};
                paramLabel = self.params{paramIdx,2};
                paramVals  = squeeze(self.sweepData(paramIdx,:,:));

                sweepPos = self.sweepRes{1}.options.sweepPosRight;

                % Plot sweep values with error bars
                errorbar(sweepPos, mean(paramVals), std(paramVals), '.b');

                % Highlight the cut-off point
                [~, cutoffIdx] = self.findCutOffForce();
                hold('on');
                plot(sweepPos(cutoffIdx), mean(paramVals(:,cutoffIdx)), '.r', 'MarkerSize', 16);

                xlim(self.sweepBoundaries);
                if ~isempty(args.ylims)
                    ylim(args.ylims(paramIdx,:));
                end

                xlabel('Fmax (pN)');
                ylabel(paramLabel);

                title(paramName);
            end
            % << nested function

            % >> nested function
            function makeErrPlot()
                paramVals = zeros(self.nBootstrapIter, self.nSweeps);

                for i = 1:self.nBootstrapIter
                    for j = 1:self.nSweeps
                        paramVals(i,j) = self.sweepRes{i}.params.err(j);
                    end
                end

                sweepPos = self.sweepRes{1}.options.sweepPosRight;

                plot(sweepPos, mean(paramVals), '.b');
                [~, cutoffIdx] = self.findCutOffForce();
                hold('on');
                plot(sweepPos(cutoffIdx), mean(paramVals(:,cutoffIdx)), '.r', 'MarkerSize', 16);

                xlim(self.sweepBoundaries);

                xlabel('Fmax (pN)');
                ylabel('err');

                title('err');
            end
            % << nested function
        end

    end % methods

    methods (Access = private)

        function [sweepResults, sweepData] = analyzeData(self, data)
            % ANALYZEDATA Function that does the actual analysis work.

            if isempty(self.sweepBoundaries)
                disp('Auto-setting sweep boundaries...');
                self.sweepBoundaries = [40 max(data.concatenatedData().f)];
            end

            % Run N bootstrap iterations, for each iteration doing an "Fmax"
            % sweep on a resampled dataset.
            disp('Bootstrapping & sweeping...');
            iterRes = cell(self.nBootstrapIter,1);
            nDatasets = data.length;

            hWait = waitbar(0, 'Bootstrapping & sweeping...');
            hWait_cleanup = onCleanup(@() close(hWait));

            for i = 1:self.nBootstrapIter
                waitbar(i/self.nBootstrapIter, hWait);

                % Get resampled dataset (or use original data for 1st iteration).
                if i == 1
                    curData = data;
                else
                    curData = getResampledDataset(data);
                end

                curDataCombined = curData.concatenatedData();

                iterRes{i} = sweepfdfits(...
                                      curDataCombined ...
                                    , 'model',              'twlc-lc-fixed' ...
                                    , 'sweepType',          'right' ...
                                    , 'sweepBoundaryType',  'f' ...
                                    , 'nSweeps',            self.nSweeps ...
                                    , 'sweepBoundaries',    self.sweepBoundaries ...
                                    , 'fitfdArgs',          {'lBounds',      self.twlcFitBounds(1,:), ...
                                                             'uBounds',      self.twlcFitBounds(2,:), ...
                                                             'startParams',  self.twlcFitStartParams, ...
                                                             'C',            self.C, ...
                                                             'Fc',           self.Fc } ...
                                    , 'Lc',                 self.Lc ...
                                    );
            end

            sweepResults = iterRes;
            sweepData    = convertSweepRes(iterRes);

            % >> nested function
            function [resampledData] = getResampledDataset(data)
                % Returns an FdDataCollection with a resampled bootstrap
                % dataset, created by resampling from the original (aligned)
                % data, with replacement.
                resampledData = FdDataCollection();
                for k = datasample(1:nDatasets, nDatasets)
                    resampledData.add('skipduplicatecheck', data.items{k});
                end
            end
            % << nested function

            % >> nested function
            function [convSweepRes] = convertSweepRes(sweepRes)
                sweepPos = sweepRes{1}.options.sweepPosRight;
                nParams  = length(sweepRes{1}.params.all{1});

                convSweepRes = zeros(nParams, self.nBootstrapIter, length(sweepPos));
                    % sweep results, re-organized as a 3D matrix:
                    %   convSweepRes(parameter, bootstrap-iter, sweep)
                for iBootstrapIter = 1:self.nBootstrapIter
                    for iSweep = 1:length(sweepPos)
                        convSweepRes(:,iBootstrapIter,iSweep) = ...
                            sweepRes{iBootstrapIter}.params.all{iSweep};
                    end
                end
            end
            % << nested function
        end

    end % methods (Access = private)

end
