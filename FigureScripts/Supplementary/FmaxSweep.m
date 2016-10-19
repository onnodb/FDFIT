classdef FmaxSweep < Figure
    % Showing the fit results for a particular magnesium concentration as a
    % function of the maximum force F[max].

    properties

        % TwlcVsMgAnalysis object. Can be left empty, in which case the
        % object will be constructed as necessary. If you already have
        % a TwlcVsMgAnalysis object containing analysis results, you can
        % also pass it in during construction of the figure:
        %
        %   >> f = Figure1('tvma', myTwlcVsMgAnalysis);
        tvma = [];

        % For which magnesium concentration to make the plot (in mM).
        mgConc = 0;  % mM

    end

    methods

        function [self] = FmaxSweep(varargin)
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

            mgConcIdx = self.tvma.findMgConcIdx(self.mgConc);

            curSweepRes  = self.tvma.sweepRes{mgConcIdx};
            curSweepData = self.tvma.sweepData{mgConcIdx};

            %% Create plot
            nParams = size(self.tvma.params,1) + 1;
            nRows = floor(sqrt(nParams));
            nCols = ceil(nParams/nRows);

            for iParam = 1:nParams
                subplot(nRows,nCols,iParam);
                self.setupAxes();
                if iParam < nParams
                    makeParamPlot(iParam);
                else
                    makeErrPlot();
                end
            end

            % >> nested functions
            function makeParamPlot(paramIdx)
                paramName  = self.tvma.params{paramIdx,1};
                paramLabel = self.tvma.params{paramIdx,2};
                paramVals  = squeeze(curSweepData(paramIdx,:,:));

                sweepPos = curSweepRes{1}.options.sweepPosRight;

                % Plot sweep values with error bars
                errorbar(sweepPos, mean(paramVals), std(paramVals), '.b');

                % Highlight the cut-off point
                [~, cutoffIdx] = self.tvma.findCutOffForce(self.mgConc);
                hold('on');
                plot(sweepPos(cutoffIdx), mean(paramVals(:,cutoffIdx)), ...
                    'dr', 'MarkerSize', 6, 'MarkerFaceColor', [1 0 0]);

                xlim(self.tvma.sweepBoundaries);

                xlabel('F_{max} (pN)');
                ylabel(paramLabel);

                title(paramName);
            end
            function makeErrPlot()
                paramVals = zeros(self.tvma.nBootstrapIter, self.tvma.nSweeps);

                for i = 1:self.tvma.nBootstrapIter
                    for j = 1:self.tvma.nSweeps
                        paramVals(i,j) = curSweepRes{i}.params.err(j);
                    end
                end

                sweepPos = curSweepRes{1}.options.sweepPosRight;

                plot(sweepPos, mean(paramVals), '.b');
                [~, cutoffIdx] = self.tvma.findCutOffForce(self.mgConc);
                hold('on');
                plot(sweepPos(cutoffIdx), mean(paramVals(:,cutoffIdx)), '.r', 'MarkerSize', 16);

                xlim(self.tvma.sweepBoundaries);

                xlabel('F_{max} (pN)');
                ylabel('R^2');

                title('R^2');
            end
            % << nested functions
        end

    end

end
