classdef GvsMagnesium < Figure
    % Graphs of g(F) for all magnesium concentrations.

    properties

        % TwlcVsMgAnalysis object. Can be left empty, in which case the
        % object will be constructed as necessary. If you already have
        % a TwlcVsMgAnalysis object containing analysis results, you can
        % also pass it in during construction of the figure:
        %
        %   >> f = Figure1('tvma', myTwlcVsMgAnalysis);
        tvma = [];

        % Change to set plot colors manually (Mx3 matrix with R,G,B colors).
        plotColors         = [];

    end

    methods

        function [self] = GvsMagnesium(varargin)
            self = self@Figure(varargin{:});
        end

        function analyze(self)
            % ANALYZE Perform calculations.

            %% Remove any data that doesn't have proper overstretching data
            % (we can't align it, and thus not analyze it)
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
        end

        function plot(self, varargin)
            % PLOT Create plot.

            defArgs = struct(...
                              'fRange',             [0 70] ...
                            , 'Fc',                 30.6 ...
                            );
            args = parseArgs(varargin, defArgs);

            plot@Figure(self);

            %% Set color defaults, if no plot colors were given.
            if isempty(self.plotColors)
                self.plotColors = copper(length(self.tvma.mgConcs));
            end

            %% Get data to plot.
            plotVals = self.tvma.getMgVsTwlcPlotData();

            Fvals = linspace(args.fRange(1), args.fRange(2), 500);
            legends = cell(size(self.tvma.mgConcs));

            for iMgConc = 1:length(self.tvma.mgConcs)
                curMg = self.tvma.mgConcs(iMgConc);

                g0 = plotVals(3,iMgConc,1);
                g1 = plotVals(4,iMgConc,1);

                plot(Fvals, n_gFun(Fvals, g0, g1, args.Fc), '-', ...
                     'Color', self.plotColors(iMgConc,:), 'LineWidth', 2);
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

    end

end
