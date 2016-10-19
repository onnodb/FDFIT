classdef MgVsOSPlateauForce < Figure
    % Average overstretching plateau force, vs. magnesium concentration.

    properties

        % Relevant magnesium concentrations (mM).
        mgConcs = [0 25 50 70 80 100 150];

        % Which distance region to use for calculating the OS plateau force.
        % Vector [min max] (um).
        OSPlateauRegion = [18.5 20];

    end

    properties (SetAccess=private)

        % Calculated overstretching plateau forces.
        OSPlateauForces = [];

    end

    methods

        function [self] = MgVsOSPlateauForce(varargin)
            self = self@Figure(varargin{:});
        end

        function analyze(self)
            % ANALYZE Perform calculations.

            %% Remove any data that doesn't have proper overstretching data.
            self.data.remove(self.data.getByTag('!nooverstretchingdata'));

            %% Calculate average OS plateau forces per concentration.
            self.OSPlateauForces = zeros(length(self.mgConcs),2);
                % col 1 = mean; col 2 = SEM

            for i = 1:length(self.mgConcs)
                curMgConc = self.mgConcs(i);

                curMgConcData = self.data.getByTag(mgConcToFdTag(curMgConc));

                forcesPerMol = zeros(curMgConcData.length,1);
                for j = 1:curMgConcData.length
                    fd = curMgConcData.items{j};
                    forcesPerMol(j) = mean(fd.f(fd.d >= self.OSPlateauRegion(1) & ...
                                                fd.d <= self.OSPlateauRegion(2)));
                end

                self.OSPlateauForces(i,1) = mean(forcesPerMol);
                self.OSPlateauForces(i,2) = std(forcesPerMol)/sqrt(curMgConcData.length);
            end

        end

        function plot(self)
            % PLOT Create plot.

            plot@Figure(self);

            self.setupAxes();

            errorbar(self.mgConcs, self.OSPlateauForces(:,1), self.OSPlateauForces(:,2), '.b', 'MarkerSize', 20);

            xlim([min(self.mgConcs)-5 max(self.mgConcs)+5]);

            set(gca(), 'YTick', floor(min(self.OSPlateauForces(:,1))):1:ceil(max(self.OSPlateauForces(:,1))));

            xlabel('[Mg] (mM)');
            ylabel('F_{os} (pN)');
        end

    end

end
