classdef Figure3 < Figure
    % Zoom-in of the transition into B-to-S overstretching for
    % different magnesium concentrations.

    properties

        % ----- Parameters

        % Relevant magnesium concentrations.
        % Vector with concentrations, in mM.
        mgConcs            = [0 25 50 100];

        % Change to set plot colors manually (Mx3 matrix with R,G,B colors).
        plotColors         = [];

        % Optional additional arguments for the 'averageFdData' function.
        averageFdData_args = {};

        % Distance region for the overstretching plateau.
        % Should be a vector [Dmin Dmax], with Dmin and Dmax both
        % in um.
        osPlateauRegion                         = [18.5 20];


        % ----- Cache

        % Aligned versions of the input data. Will be generated
        % automatically if not given. If you have already performed
        % alignment yourself, you can also pass in this aligned version
        % by using this property.
        % If you have a "TwlcVsMgAnalysis" object with analysis results
        % available, you can use the "dataAligned" property of the
        % "TwlcVsMgAnalysis" object, to save yourself from having to do
        % duplicate work.
        dataAligned                             = [];

    end

    properties (SetAccess = private)

        % Averaged data curves, per Mg concentration.
        dataAveraged       = [];

        % Axes handle to the main plot.
        axes_outer         = [];

    end

    methods

        function [self] = Figure3(varargin)
            self = self@Figure(varargin{:});
        end

        function analyze(self)
            % ANALYZE Perform calculations.

            %% Remove any data that doesn't have proper overstretching data
            % (which is necessary for vertical alignment)
            self.data.remove(self.data.getByTag('!nooverstretchingdata'));

            %% Align the F,d curves (all together).
            if isempty(self.dataAligned)
                disp('Aligning...');
                self.dataAligned = alignFdCurves(self.data, 'OSPlateauRegion', self.osPlateauRegion);
            else
                disp('Using pre-set aligned data.');
            end

            %% Create averaged datasets.
            self.dataAveraged = FdDataCollection();
            for i = 1:length(self.mgConcs)
                curMgConc = self.mgConcs(i);

                self.dataAveraged.add(...
                        averageFdData(self.dataAligned.getByTag(mgConcToFdTag(curMgConc)), ...
                        self.averageFdData_args{:}) ...
                        );
            end
        end

        function plot(self, plotMgConcs)
            % PLOT Create plot.

            plot@Figure(self);

            if nargin < 2
                plotMgConcs = self.mgConcs;
            end

            %% Set color defaults, if no plot colors were given.
            if isempty(self.plotColors)
                self.plotColors = hsv(length(plotMgConcs));
            end

            %% Create main plot
            self.axes_outer = axes();
            self.setupAxes();

            legendItems = cell(length(plotMgConcs),1);
            for i = 1:length(plotMgConcs)
                itemIdx = find(self.mgConcs==plotMgConcs(i), 1, 'first');
                assert(~isempty(itemIdx));

                curFd = self.dataAveraged.items{itemIdx};

                plot(self.axes_outer, curFd.d, curFd.f, '-', 'Color', self.plotColors(i,:), 'LineWidth', 2);

                legendItems{i} = sprintf('%d mM', self.mgConcs(itemIdx));
                if i == 1
                    hold(self.axes_outer, 'on');
                end
            end

            xlim(self.axes_outer, [16.3 18]);
            ylim(self.axes_outer, [30 70]);

            set(self.axes_outer, 'XTick', 16:0.5:18);
            set(self.axes_outer, 'YTick', 30:10:70);

            xlabel(self.axes_outer, 'Distance ({\mu}m)');
            ylabel(self.axes_outer, 'Force (pN)');

            legend(self.axes_outer, legendItems{:});
            legend(self.axes_outer, 'Location', 'NorthWest');
        end

    end % methods

end
