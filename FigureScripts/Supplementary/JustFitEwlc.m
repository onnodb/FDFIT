classdef JustFitEwlc < Figure
    % Fit only the eWLC to the data, using different approaches, and compare
    % them to the results from the full tWLC analysis.

    properties

        % TwlcVsMgAnalysis object. Can be left empty, in which case the
        % object will be constructed as necessary. If you already have
        % a TwlcVsMgAnalysis object containing analysis results, you can
        % also pass it in during construction of the figure:
        %
        %   >> f = LpAndSAlternatives('tvma', myTwlcVsMgAnalysis);
        tvma                    = [];

        % Contour length of the DNA (in um).
        Lc = 16.5;

        % Default, global plot arguments.
        globalPlotArgs = {'MarkerSize', 20};

        % Extra plot settings for each of the items.
        plotArgs = {  {}                                                                             ... % 1. raw data, individual fits
                    , {'Marker', 'o', 'MarkerSize', 6}                                               ... % 2. fit combined aligned data
                    , {'Marker', 'd', 'MarkerSize', 6, 'MarkerFaceColor', 'auto', 'LineStyle', '--'} ... % 3. tWLC analysis
                    };

        % How many bootstrap iterations do perform, for getting error bars
        % on global fits?
        nBootstrapIter = 100;

    end % properties

    properties (SetAccess=private)

        items = {...
                  'Raw data, individual fits'               ...
                ; 'Fit combined aligned data'               ...
                ; 'tWLC analysis'                           ...
                };

        %% --- Calculated values

        LpSData;

    end % properties

    % ------------------------------------------------------------------------

    methods

        function [self] = JustFitEwlc(varargin)
            self = self@Figure(varargin{:});
        end

        function analyze(self)
            % ANALYZE Perform calculations.

            %% Run / load TwlcVsMgAnalysis.
            if isempty(self.tvma)
                self.tvma = TwlcVsMgAnalysis('data', self.data);
                self.tvma.analyze();
            elseif ischar(self.tvma)
                fprintf('Loading TwlcVsMgAnalysis object from: %s\n', self.tvma);
                self.tvma = load(self.tvma);
            elseif isa(self.tvma, 'TwlcVsMgAnalysis')
                disp('Using pre-set TwlcVsMgAnalysis object.');
            else
                error('Invalid value for "tvma" property.');
            end

            self.calculate();
        end

        function calculate(self, whichItems, forceRecalculate)
            % CALCULATE Calculate values for different analysis approaches ("items").

            if nargin < 2
                whichItems = 1:length(self.items);
            end
            if nargin < 3
                forceRecalculate = false;
            end

            if isempty(self.LpSData)
                nItems = length(self.items);
                self.LpSData = struct('Lp', {}, 'S', {});
                self.LpSData(nItems).Lp = [];
                for i = 1:nItems
                    self.LpSData(i).Lp = [];
                    self.LpSData(i).S  = [];
                end
            end

            assert(~isempty(whichItems) && isvector(whichItems) && ...
                   all(whichItems>0 & whichItems<=length(self.items)));

            for itemIdx = whichItems(:)'
                if ( isempty(self.LpSData(itemIdx).Lp) || isempty(self.LpSData(itemIdx).S) ) ...
                        || forceRecalculate
                    fprintf('\n>>>>>>>>>> Calculating item %d: %s...\n', itemIdx, self.items{itemIdx});
                    switch itemIdx
                        case 1
                            %% -- 1. Raw data, individual fits
                            [Lp, S] = self.getFromIndividualFits(self.tvma.data);
                        case 2
                            %% -- 2. Simple fits to combined aligned data
                            [Lp, S] = self.getFromCombinedAlignedDataFits();
                        case 3
                            %% -- 3. tWLC analysis results
                            twlcData = self.tvma.getMgVsTwlcPlotData();
                            Lp = squeeze(twlcData(1,:,:));
                            S  = squeeze(twlcData(2,:,:));
                        otherwise
                            error('Internal error: invalid item index %d.', itemIdx);
                    end
                    self.LpSData(itemIdx).Lp = Lp;
                    self.LpSData(itemIdx).S  = S;
                end
            end % for itemIdx
        end

        function plot(self, varargin)
            % PLOT Create plot.

            defArgs = struct(...
                              'items',                  [] ...
                            , 'plotArgs',               {{}} ...
                            , 'stagger',                true ...
                            , 'staggerOffset',          1.5 ...
                            );
            args = parseArgs(varargin, defArgs, {'stagger'});

            plot@Figure(self);

            if isempty(args.items)
                args.items = 1:length(self.items);
            end

            subplot(1,2,1);
            self.setupAxes();
            n_makePlot('Lp');

            subplot(1,2,2);
            self.setupAxes();
            n_makePlot('S');

            legends = self.items(args.items);
            legend(legends{:});

            % >> nested functions
            function n_makePlot(param)
                staggerIdx = 1;
                for itemIdx = args.items(:)'
                    paramVals = self.LpSData(itemIdx).(param);

                    curPlotArgs = [self.globalPlotArgs self.plotArgs{itemIdx} args.plotArgs];

                    if args.stagger
                        mgConcs = self.tvma.mgConcs + ((staggerIdx-length(args.items)/2)*args.staggerOffset);
                    else
                        mgConcs = self.tvma.mgConcs;
                    end

                    if size(paramVals,2) == 1
                        plot(mgConcs, paramVals(:,1), '.', curPlotArgs{:});
                    elseif size(paramVals,2) == 2
                        errorbar(mgConcs, paramVals(:,1), paramVals(:,2), '.', curPlotArgs{:});
                    elseif size(paramVals,2) == 3
                        errorbar(mgConcs, paramVals(:,1), paramVals(:,2), paramVals(:,3), '.', curPlotArgs{:});
                    end
                    hold('on');

                    staggerIdx = staggerIdx + 1;
                end
                xlim([min(self.tvma.mgConcs)-10 max(self.tvma.mgConcs)+10]);
                xlabel('[Mg] (mM)');
                ylabel(param);
            end
            % << nested functions
        end

    end

    methods

        function [Lp, S] = getFromIndividualFits(self, data, varargin)
            defArgs = struct(...
                              'makePlots',                  false ...
                            , 'inspectFits',                false ...
                            );
            args = parseArgs(varargin, defArgs, {'makePlots', 'inspectFits'});

            Lp = zeros(length(self.tvma.mgConcs),2);
            S  = zeros(length(self.tvma.mgConcs),2);
            for iMgConc = 1:length(self.tvma.mgConcs)
                curMgConc = self.tvma.mgConcs(iMgConc);

                mgData = data.getByTag(mgConcToFdTag(curMgConc));
                assert(~mgData.isempty());
                mgFits = fitfd(mgData, 'model', 'odijk-d0-f0', 'Lc', self.Lc);

                Lp(iMgConc,:) = [mean([mgFits.Lp]) std([mgFits.Lp])/sqrt(length(mgFits))];  % mean +/- SEM
                S(iMgConc,:)  = [mean([mgFits.S ]) std([mgFits.S ])/sqrt(length(mgFits))];  % mean +/- SEM

                if args.makePlots
                    n_makePlot(mgFits, Lp(iMgConc), S(iMgConc), curMgConc);
                end
                if args.inspectFits
                    explorefdfits(mgFits);
                end
            end

            % >> nested functions
            function n_makePlot(mgFits, Lp, S, mgConc)
                figure();
                subplot(1,2,1);
                    plotkde([mgFits.Lp], 'normalize');
                    hold('on');
                    plot(Lp, 1, '.r', 'MarkerSize', 20);
                    xlabel('Lp');
                subplot(1,2,2);
                    plotkde([mgFits.S], 'normalize');
                    hold('on');
                    plot(S, 1, '.r', 'MarkerSize', 20);
                    xlabel('S');

                title(sprintf('[Mg] = %d mM', mgConc));
            end
            % << nested functions
        end

        function [Lp, S] = getFromCombinedAlignedDataFits(self)
            Lc_ = self.Lc;

            Lp = zeros(length(self.tvma.mgConcs),2);
            S  = zeros(length(self.tvma.mgConcs),2);

            for iMgConc = 1:length(self.tvma.mgConcs)
                curMgConc = self.tvma.mgConcs(iMgConc);
                fprintf('\n##### [Mg] = %d mM\n Bootstrapping...\n', curMgConc);

                mgData = self.tvma.dataAligned.getByTag(mgConcToFdTag(curMgConc));
                assert(~mgData.isempty());

                bsLpVals = zeros(self.nBootstrapIter,1);
                bsSVals  = zeros(self.nBootstrapIter,1);
                parfor iBootstrapIter = 1:self.nBootstrapIter
                    if iBootstrapIter == 1
                        bsData = mgData;
                    else
                        bsData = resampleFdDataCollection(mgData);
                    end
                    f = fitfd(bsData.concatenatedData(), 'model', 'odijk-lc-fixed', 'Lc', Lc_);
                    bsLpVals(iBootstrapIter) = f.Lp;
                    bsSVals(iBootstrapIter)  = f.S;
                end

                if self.nBootstrapIter > 1
                    Lp(iMgConc,:) = [mean(bsLpVals) std(bsLpVals)];
                    S(iMgConc,:)  = [mean(bsSVals)  std(bsSVals) ];
                else
                    Lp(iMgConc,:) = [bsLpVals NaN];
                    S(iMgConc,:)  = [bsSVals  NaN];
                end
            end
        end

    end % methods

end
