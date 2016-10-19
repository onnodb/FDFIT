classdef Figure1 < Figure
    % Comparing an eWLC fit to a tWLC fit.
    % At the same time, shows one of the final fits to an averaged data curve.

    properties

        % TwlcVsMgAnalysis object. Can be left empty, in which case the
        % object will be constructed as necessary. If you already have
        % a TwlcVsMgAnalysis object containing analysis results, you can
        % also pass it in during construction of the figure:
        %
        %   >> f = Figure1('tvma', myTwlcVsMgAnalysis);
        tvma                    = [];

        % Use data at this magnesium concentration.
        mgConc                  = 0;

        % ----- Optional

        % FdData object for additional plot in inset.
        % If not set, automatically tries to find an FdData object
        % tagged "lowsalt" in the TwlcVsMgAnalysis's "data" property.
        lowSaltCurve            = [];

        % Cell array with any additional arguments to pass into the
        % "averageFdData" function (during "analyze").
        averageFdData_args      = {};

    end

    properties (SetAccess = private)

        % Holds data that has been averaged during the call to "analyze".
        dataAveraged            = [];

        % Fit results of the eWLC model.
        fitEWLC                 = [];

        % Fit results of the tWLC model.
        fitTWLC                 = [];

        % Axes handle for main plot.
        axes_outer              = [];

        % Axes handle for inset plot.
        axes_inset              = [];

    end % properties

    % ------------------------------------------------------------------------

    methods

        function [self] = Figure1(varargin)
            self = self@Figure(varargin{:});
        end

        function analyze(self)
            % ANALYZE Perform calculations.

            %% Remove any data that doesn't have proper overstretching plateaus.
            % (We can't align it, and thus not analyze it).
            self.data.remove(self.data.getByTag('!nooverstretchingdata'));

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

            %% Get averaged data curve.
            disp('Calculating averaged data curve...');
            self.dataAveraged = averageFdData(...
                                    self.tvma.dataAligned.getByTag(mgConcToFdTag(self.mgConc)), ...
                                    self.averageFdData_args{:} ...
                                    );

            %% Fit eWLC.
            disp('Fitting eWLC...');
            fitE = fitfd(self.dataAveraged, 'model', 'odijk');
            self.fitEWLC = FdData('eWLC fit', fitE(self.dataAveraged.d), self.dataAveraged.d, 1:length(self.dataAveraged.d));
            disp(fitE);

            %% Fit curve for tWLC.
            disp('Calculating tWLC fit curve...');
            iMgConc = find(self.tvma.mgConcs==self.mgConc, 1, 'first');
            assert(~isempty(iMgConc));
            tvmaVals = self.tvma.getMgVsTwlcPlotData();  % (iParam, iMgConc, 1=mean/2=errorbar)
            self.fitTWLC = FdData('tWLC fit', ...
                                  self.dataAveraged.f, ...
                                  ftWLC(...
                                        self.dataAveraged.f, ...
                                        tvmaVals(1,iMgConc,1), 16.5, tvmaVals(2,iMgConc,1), ...
                                        440, tvmaVals(3,iMgConc,1), tvmaVals(4,iMgConc,1), 30.6 ...
                                        ), ...
                                  (1:length(self.dataAveraged.f)) ...
                                  );

            %% Get extra low-salt curve for inset.
            if isempty(self.lowSaltCurve)
                disp('Fetching low-salt curve...');
                lowSalt = self.data.getByTag('lowsalt');
                if ~lowSalt.isempty()
                    self.lowSaltCurve = lowSalt.items{1};
                end
            end

            %% Done!
            disp('Done.');
        end

        function plot(self)
            % PLOT Create plot.

            plot@Figure(self);

            %% Main plot: overall plot
            self.axes_outer = axes();
            self.setupAxes();
            plotDataAndFits(self.axes_outer, true);

            % Tweak axes appearance
            xlabel(self.axes_outer, 'Distance ({\mu}m)');
            ylabel(self.axes_outer, 'Force (pN)');
            xlim(self.axes_outer, [13 19]);
            ylim(self.axes_outer, [0 70]);
            set(self.axes_outer, 'XTick', 13:2:19);
            set(self.axes_outer, 'YTick', 0:20:70);


            %% Inset: zoom in of difference between eWLC / tWLC
            self.axes_inset = axes('Position', [0.21 0.52 0.32 0.35]);
            self.setupAxes();
            plotDataAndFits(self.axes_inset, false);

            % Tweak axes appearance
            xlim(self.axes_inset, [16.3 17.6]);
            ylim(self.axes_inset, [30 70]);
            set(self.axes_inset, 'XTick', [16.5 17.5]);
            set(self.axes_inset, 'YTick', [30 50 70]);


            function plotDataAndFits(ax, plotLowSaltCurve)
                % Extra inset curve, if requested
                if ~isempty(self.lowSaltCurve) && plotLowSaltCurve
                    plot(ax, self.lowSaltCurve.d, self.lowSaltCurve.f, '-.', 'Color', [0.6 0.6 0.6]);
                    hold(ax, 'on');
                end

                % Data curve
                plot(ax, self.dataAveraged.d(1:3:end), self.dataAveraged.f(1:3:end), 'ok', 'MarkerSize', 5);

                hold(ax, 'on');

                % eWLC fit
                hEWLCPlot = plot(ax, self.fitEWLC.d, self.fitEWLC.f, '--r', 'LineWidth', 2);

                % tWLC fit
                hTWLCPlot = plot(ax, self.fitTWLC.d, self.fitTWLC.f, '-', 'Color', [0.3 0.3 1], 'LineWidth', 2);
            end

        end

    end % methods

end
