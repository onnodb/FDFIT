classdef Figure5 < Figure
    % tWLC parameters as a function of Magnesium concentration.

    properties

        % TwlcVsMgAnalysis object. Can be left empty, in which case the
        % object will be constructed as necessary. If you already have
        % a TwlcVsMgAnalysis object containing analysis results, you can
        % also pass it in during construction of the figure:
        %
        %   >> f = Figure5('tvma', myTwlcVsMgAnalysis);
        tvma = [];

    end % properties

    % ------------------------------------------------------------------------

    methods

        function [self] = Figure5(varargin)
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

        function plot(self)
            % PLOT Create plot.

            plot@Figure(self);

            plotVals = self.tvma.getMgVsTwlcPlotData();

            nParams = size(self.tvma.params,1);
            nRows = floor(sqrt(nParams));
            nCols = ceil(nParams/nRows);

            for iParam = 1:nParams
                subplot(nRows,nCols,iParam);
                self.setupAxes();

                errorbar(self.tvma.mgConcs, plotVals(iParam,:,1), plotVals(iParam,:,2), '.b', 'MarkerSize', 16);

                xlim([min(self.tvma.mgConcs)-5 max(self.tvma.mgConcs)+5]);

                maxAbsVal = max(abs(plotVals(iParam,:,1))) * sign(plotVals(iParam,1,1));
                minAbsVal = min(abs(plotVals(iParam,:,1))) * sign(plotVals(iParam,1,1));
                ylim(sort([minAbsVal/1.05 maxAbsVal*1.05]));

                xlabel('[Mg] (mM)');
                ylabel(self.tvma.params{iParam,2});
            end

        end

    end % methods

end
