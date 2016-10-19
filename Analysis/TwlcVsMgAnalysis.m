classdef TwlcVsMgAnalysis < handle
    % TWLCVSMGANALYSIS Encapsulation of the analysis of the tWLC for different Mg concentrations.
    %
    % This is a convenience class that encapsulates the analysis procedure
    % for the paper: extracting the dependence of tWLC fit parameters
    % on magnesium concentration.
    %
    % The idea: set a bunch of properties to pass in the data and algorithm
    % parameters, run "analyze", and store the resulting object on disk
    % for later retrieval. This object can now be passed into the various
    % "FigureN" classes to re-create the figures from the paper.

    properties

        % ----- Input data

        % FdDataCollection containing the data to be analyzed.
        % Each FdData object should contain a clean F,d curve, and should
        % be tagged with a tag following the format of the "fdTagToMgConc"
        % and "mgConcToFdTag" functions (see under "Utils/"). The tags are
        % used to link the F,d data to magnesium concentrations.
        data                                    = [];


        % ----- Cache

        % Aligned versions of the input data. Will be generated
        % automatically if not given. If you have already performed
        % alignment yourself, you can also pass in this aligned version
        % by using this property.
        dataAligned                             = [];


        % ----- Parameters

        % Magnesium concentrations to analyze.
        % Vector of concentrations, in mM.
        mgConcs                                 = [0 25 50 70 80 100 150];


        % -- Sweeping & bootstrapping

        % Contour length of the DNA (um).
        Lc                                      = 16.5;

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
        twlcFitBounds                           = [5 500 -2000 2; ... % Lp S g0 g1
                                                   150 5000 -100 60];

        % Optional vector with starting values for the fit parameters
        % [Lp S g0 g1].
        twlcFitStartParams                      = [];


        % -- Misc.

        % Parameter names (for plotting; usually no need to change this).
        params = {  'Lp',       'L_p (nm)',     'nm'        , [41 43]     ; ...
                    'S',        'S (pN)',       'pN'        , []          ; ...
                    'g0',       'g_0 (pN nm)',  'pN nm'     , [-600 -350] ; ...
                    'g1',       'g_1 (nm)',     'nm'        , [12 17]       ...
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

        function [self] = TwlcVsMgAnalysis(varargin)
            % TWLCVSMGANALYSIS Constructor.
            %
            % SYNTAX:
            % tvma = TwlcVsMgAnalysis('key', value, ...)
            %
            % KEY-VALUE PAIR ARGUMENTS:
            % All class properties can be initialized using the key-value pair
            % system. (Note: these are case-sensitive).

            parseClassArgs(varargin, self);

            if isempty(self.data) || ~isa(self.data, 'FdDataCollection') || (self.data.length == 0)
                error('No data or invalid data given: FdDataCollection expected');
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


            %% Analyze data per concentration.
            nMgConcs = length(self.mgConcs);

            self.sweepRes    = cell(nMgConcs,1);
            self.sweepData   = cell(nMgConcs,1);

            for i = 1:length(self.mgConcs)
                curMgConc = self.mgConcs(i);
                fprintf('\n##### [Mg] = %d mM\n\n', curMgConc);

                mgData = self.dataAligned.getByTag(mgConcToFdTag(curMgConc));

                if mgData.isempty()
                    error('No data for magnesium concentration %d mM.', curMgConc);
                end

                [self.sweepRes{i}, self.sweepData{i}] = self.analyzeData(mgData);
            end
        end

        function [cutoffForce, cutoffIdx] = findCutOffForce(self, mgConc)
            % FINDCUTOFFFORCE Finds the cut-off force for a particular Mg concentration.
            %
            % Determines the upper limit of validity of the tWLC, for a
            % particular Mg concentration.
            %
            % INPUT:
            % mgConc = Mg concentration (mM); must be present in the list
            %       'self.mgConcs'.
            %
            % OUTPUT:
            % cutoffForce = upper limit of validity of the tWLC, i.e., the
            %       maximum force that can be used for fitting the tWLC (pN).
            % cutoffIdx = index of the sweep that corresponds to 'cutoffForce'
            %       (index into 'sweepRes{mgIdx}.options.sweepPosRight).

            if isempty(self.sweepData)
                error('No analysis data; run "analyze" first.');
            end

            mgConcIdx = self.findMgConcIdx(mgConc);

            curSweepRes  = self.sweepRes{mgConcIdx};
            curSweepData = self.sweepData{mgConcIdx};

            % Cut-off force determined from maximum GOF (goodness-of-fit).
            gofVals = zeros(self.nBootstrapIter, self.nSweeps);
            for i = 1:self.nBootstrapIter
                for j = 1:self.nSweeps
                    gofVals(i,j) = curSweepRes{i}.params.err(j);
                end
            end
            gofVals = mean(gofVals);
            [~, maxGofIdx] = max(gofVals);
            cutoffForce = curSweepRes{1}.options.sweepPosRight(maxGofIdx);
            cutoffIdx   = maxGofIdx;
        end

        function [plotVals] = getMgVsTwlcPlotData(self)
            % GETMGVSTWLCPLOTDATA Collects plot data for 'plotMgVsTwlc'.

            nParams  = size(self.params,1);
            nMgConcs = length(self.mgConcs);

            plotVals = zeros(nParams, nMgConcs, 2);
                % plot values: (iParam, iMgConc, 1=mean/2=errorbar)

            for iMgConc = 1:nMgConcs
                curMgConc       = self.mgConcs(iMgConc);
                [~, curCutoff]  = self.findCutOffForce(curMgConc);
                curSweepData    = self.sweepData{iMgConc};

                for iParam = 1:nParams
                    paramVals = squeeze(curSweepData(iParam,:,curCutoff));
                    plotVals(iParam, iMgConc, 1) = mean(paramVals);
                    plotVals(iParam, iMgConc, 2) = std(paramVals);
                end
            end
        end

        function plotMgVsCutoffForce(self)
            % PLOTMGVSCUTOFFFORCE Creates a plot of [Mg] vs Fmax (cut-off force).
            cutoffForces = zeros(size(self.mgConcs));
            for i = 1:length(self.mgConcs)
                cutoffForces(i) = self.findCutOffForce(self.mgConcs(i));
            end

            figure;
            plot(self.mgConcs, cutoffForces, '.b', 'MarkerSize', 20);
            xlabel('[Mg] (mM)');
            ylabel('Cut-off force (pN)');
        end

        function plotMgVsTwlc(self, varargin)
            % PLOTMGVSTWLC Plots the final [Mg] vs tWLC graphs.
            %
            % KEY-VALUE PAIR ARGUMENTS:
            % hFig = optional Figure handle of the figure window that should
            %       receive the plot.

            defArgs = struct(...
                              'hFig',               [] ...
                            );
            args = parseArgs(varargin, defArgs);
            if isempty(args.hFig)
                args.hFig = figure();
            else
                figure(args.hFig);
                clf();
            end

            plotVals = self.getMgVsTwlcPlotData();

            %% Create plot
            nParams = size(self.params,1);
            nRows = floor(sqrt(nParams));
            nCols = ceil(nParams/nRows);

            for iParam = 1:nParams
                subplot(nRows, nCols, iParam);

                errorbar(self.mgConcs, plotVals(iParam,:,1), plotVals(iParam,:,2), '.b');

                if ~isempty(self.params{iParam,4})
                    ylim(self.params{iParam,4});
                end

                xlabel('[Mg] (mM)');
                ylabel(self.params{iParam,2});
            end
        end

        function plotMgVsG(self, varargin)
            defArgs = struct(...
                              'hFig',               [] ...
                            , 'fRange',             [0 70] ...
                            , 'Fc',                 30.6 ...
                            );
            args = parseArgs(varargin, defArgs);

            if isempty(args.hFig)
                args.hFig = figure();
            else
                figure(args.hFig);
                clf();
            end

            plotVals = self.getMgVsTwlcPlotData();

            Fvals = linspace(args.fRange(1), args.fRange(2), 500);
            legends = cell(size(self.mgConcs));

            for iMgConc = 1:length(self.mgConcs)
                curMg = self.mgConcs(iMgConc);

                g0 = plotVals(3,iMgConc,1);
                g1 = plotVals(4,iMgConc,1);

                plot(Fvals, n_gFun(Fvals, g0, g1, args.Fc), '-');
                hold('on');
                legends{iMgConc} = sprintf('[Mg2+] = %d mM', curMg);
            end

            xlabel('F (pN)');
            ylabel('g (pN nm)');
            legend(legends{:}, 'Location', 'northwest');

            % >> nested functions
            function [g] = n_gFun(F, g0, g1, Fc)
                g = zeros(size(F));
                g(F <  Fc) = g0 + g1.*Fc;
                g(F >= Fc) = g0 + g1.*F(F>=Fc);
            end
            % << nested functions
        end

        function plotSweeps(self, mgConc, varargin)
            % PLOTSWEEPS Plot the sweep results for a particular Mg concentration.
            %
            % INPUT:
            % mgConc = Mg concentration (mM); must be present in 'self.mgConcs'
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
                            , 'ylims',                  [41 44; 1700 2000; -1400 -400; 10 35] ...
                            );
            args = parseArgs(varargin, defArgs, {'plotErr'});

            if isempty(self.sweepRes)
                error('No data to plot; run "analyze" first');
            end

            mgConcIdx = self.findMgConcIdx(mgConc);

            curSweepRes  = self.sweepRes{mgConcIdx};
            curSweepData = self.sweepData{mgConcIdx};

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

            % Add title with Mg concentration
            axes('Position', [0 0 1 1], 'Visible', 'off');
            text(0.5, 0.97, sprintf('[Mg] = %d mM', mgConc), 'Units', 'normalized', 'Parent', gca());

            % >> nested function
            function makeParamPlot(paramIdx)
                paramName  = self.params{paramIdx,1};
                paramLabel = self.params{paramIdx,2};
                paramVals  = squeeze(curSweepData(paramIdx,:,:));

                sweepPos = curSweepRes{1}.options.sweepPosRight;

                % Plot sweep values with error bars
                errorbar(sweepPos, mean(paramVals), std(paramVals), '.b');

                % Highlight the cut-off point
                [~, cutoffIdx] = self.findCutOffForce(mgConc);
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
                        paramVals(i,j) = curSweepRes{i}.params.err(j);
                    end
                end

                sweepPos = curSweepRes{1}.options.sweepPosRight;

                plot(sweepPos, mean(paramVals), '.b');
                [~, cutoffIdx] = self.findCutOffForce(mgConc);
                hold('on');
                plot(sweepPos(cutoffIdx), mean(paramVals(:,cutoffIdx)), '.r', 'MarkerSize', 16);

                xlim(self.sweepBoundaries);

                xlabel('Fmax (pN)');
                ylabel('err');

                title('err');
            end
            % << nested function
        end

        function [mgConcIdx] = findMgConcIdx(self, mgConc)
            mgConcIdx = find(self.mgConcs==mgConc);
            if isempty(mgConcIdx)
                error('No data for [Mg] = %d', mgConc);
            end
        end

    end % methods

    methods (Access = private)

        function [sweepResults, sweepData] = analyzeData(self, data)
            % ANALYZEDATA Function that does the actual analysis work.

            if isempty(self.sweepBoundaries)
                disp('Auto-setting sweep boundaries...');
                self.sweepBoundaries = [40 max(data.concatenatedData().f)];
            end

            %% Run N bootstrap iterations, for each iteration doing an "Fmax"
            % sweep on a resampled dataset
            disp('Bootstrapping & sweeping...');
            iterRes = cell(self.nBootstrapIter,1);
            nDatasets = data.length;

            hWait = waitbar(0, 'Bootstrapping & sweeping...');
            hWait_cleanup = onCleanup(@() close(hWait));

            for i = 1:self.nBootstrapIter
                waitbar(i/self.nBootstrapIter, hWait);

                % Get resampled dataset (or use original data for 1st iteration)
                if i == 1
                    curData = data;
                else
                    curData = resampleFdDataCollection(data);
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
                                                              'startParams', self.twlcFitStartParams } ...
                                    , 'Lc',                 self.Lc ...
                                    );
            end

            sweepResults = iterRes;
            sweepData    = convertSweepRes(iterRes);

            % >> nested function
            function [convSweepRes] = convertSweepRes(sweepRes)
                sweepPos = sweepRes{1}.options.sweepPosRight;
                nParams  = length(sweepRes{1}.params.all{1});

                convSweepRes = zeros(nParams, self.nBootstrapIter, length(sweepPos));
                    % sweep results, re-organized as a 3D matrix
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
