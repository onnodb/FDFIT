classdef ManualGFits < Figure
    % Manual fits of g0/g1 confirm the original analysis from the paper.

    properties

        % TwlcVsMgAnalysis object. Can be left empty, in which case the
        % object will be constructed as necessary. If you already have
        % a TwlcVsMgAnalysis object containing analysis results, you can
        % also pass it in during construction of the figure:
        %
        %   >> f = Figure1('tvma', myTwlcVsMgAnalysis);
        tvma = [];

        % Twist rigidity.
        C                       = 440;

        % Thermal energy.
        kT                      = 4.11;

        % Range of forces between which we git |g(F)| to get g0/g1.
        fitRange                = [50 60]; % pN

    end

    properties (SetAccess = private)

        % Holds data that has been averaged during the call to "analyze".
        dataCombined            = {};

        % Cell array, one cell per magnesium concentration.
        % In each cell, a struct with fields Lp and S; each is a length-2
        % vector, with the first value the fit value and the second one the
        % error bar.
        ewlcParamVals           = {};

        % |g(F)| vals, one cell per magnesium concentration; contains vector
        % with |g(F)| values.
        gVals                   = {};

        % g0/g1 fitted values.
        % (2xMx2) matrix, with indices (1=point 2=err; magnesium conc. index;
        % 1=g0, 2=g1).
        gFittedVals             = [];

    end % properties

    % ------------------------------------------------------------------------

    methods

        function [self] = ManualGFits(varargin)
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


            %% Process each of the magnesium concentrations.
            self.dataCombined = cell(size(self.tvma.mgConcs));
            self.gFittedVals = zeros(2,length(self.tvma.mgConcs),2);
            for iMgConc = 1:length(self.tvma.mgConcs)
                curMg = self.tvma.mgConcs(iMgConc);
                fprintf('Processing [Mg] = %d mM...\n', curMg);

                % Create combined curve.
                self.dataCombined{iMgConc} = self.tvma.dataAligned.getByTag(mgConcToFdTag(curMg)).concatenatedData();

                % Fit eWLC.
                f = fitfd(...
                        self.dataCombined{iMgConc}, ...
                        'model',    'odijk-lc-fixed', ...
                        'Lc',       self.tvma.Lc ...
                        );
                cf = diff(confint(f, 0.68))./2;
                self.ewlcParamVals{iMgConc} = struct(...
                          'Lp',     [f.Lp cf(1)] ...
                        , 'S',      [f.S  cf(2)] ...
                        );

                % Calculate |g(F)|.
                gSquared = ...
                        (self.ewlcParamVals{iMgConc}.S(1) * self.C) - (self.C .* self.dataCombined{iMgConc}.f) ...
                        ./ (self.dataCombined{iMgConc}.d ./ self.tvma.Lc - 1 ...
                            + 0.5 * sqrt(self.kT ./ self.dataCombined{iMgConc}.f ./ self.ewlcParamVals{iMgConc}.Lp(1)));
                gSquared(imag(sqrt(gSquared)) ~= 0) = NaN;
                curGVals = sqrt(gSquared);
                self.gVals{iMgConc} = curGVals;


                % Fit "manually".
                curFVals = self.dataCombined{iMgConc}.f;

                fitFVals = curFVals(curFVals >= self.fitRange(1) & curFVals <= self.fitRange(2));
                fitGVals = curGVals(curFVals >= self.fitRange(1) & curFVals <= self.fitRange(2));
                fitFVals(isnan(fitGVals)) = [];
                fitGVals(isnan(fitGVals)) = [];
                [f,~] = fit(fitFVals, fitGVals, 'poly1');
                ci = confint(f, 0.68);

                self.gFittedVals(:,iMgConc,1) = [f.p2 diff(ci(:,2))/2];   % g0
                self.gFittedVals(:,iMgConc,2) = [f.p1 diff(ci(:,1))/2];   % g1
            end


            %% Done!
            disp('Done.');
        end

        function plot(self)
            % PLOT Create plot.

            plot@Figure(self);

            twlcdata = self.tvma.getMgVsTwlcPlotData();

            subplot(1,3,1);
            self.plotGFit(0, 'parent', gca());

            for gIndex = 1:2
                subplot(1,3,1+gIndex);
                errorbar(...
                    self.tvma.mgConcs, self.gFittedVals(1,:,gIndex), self.gFittedVals(2,:,gIndex,:), ...
                    '.', 'MarkerSize', 24 ...
                    );
                hold('on');
                errorbar(...
                    self.tvma.mgConcs, twlcdata(2+gIndex,:,1), twlcdata(2+gIndex,:,2), ...
                    'o', 'MarkerSize', 8 ...
                    );
                xlabel('[Mg] (mM)');
                if gIndex == 1
                    ylabel('g_0 (pN nm)');
                else
                    ylabel('g_1 (nm)');
                end
                xlim([min(self.tvma.mgConcs)-10 max(self.tvma.mgConcs)+10]);
            end

            legend('Manual fit', 'Original tWLC analysis');
        end

        function plotGFit(self, mgConc, varargin)
            % PLOTGFIT Plot the fit of the |g(F)| graph, for a particular Mg concentration
            %
            % INPUT:
            % mgConc = magnesium concentration (in mM).

            defArgs = struct(...
                              'xlim',                   [0 70] ...
                            , 'ylim',                   [0 1e3] ...
                            , 'parent',                 [] ...
                            );

            args = parseArgs(varargin, defArgs);
            iMgConc = find(self.tvma.mgConcs==mgConc);
            assert(~isempty(iMgConc));

            curGVals = self.gVals{iMgConc};
            curFVals = self.dataCombined{iMgConc}.f;

            fitFVals = curFVals(curFVals >= self.fitRange(1) & curFVals <= self.fitRange(2));
            fitGVals = curGVals(curFVals >= self.fitRange(1) & curFVals <= self.fitRange(2));
            fitFVals(isnan(fitGVals)) = [];
            fitGVals(isnan(fitGVals)) = [];

            if isempty(args.parent)
                figure();
            else
                axes(args.parent);
            end
            plot(curFVals, curGVals, '.', 'Color', [0.7 0.7 0.7]);
            hold('on');
            plot(fitFVals, fitGVals, '.', 'Color', [0 0 1]);
            GFunc = @(g0,g1,F) g0 + g1.*F;
            plot(...
                fitFVals, GFunc(self.gFittedVals(1,iMgConc,1), self.gFittedVals(1,iMgConc,2), fitFVals), ...
                '-r', 'LineWidth', 1.5 ...
                );
            xlim(args.xlim);
            ylim(args.ylim);
            xlabel('F (pN)');
            ylabel('|g(F)| (pN nm)');
        end

    end % methods

end
