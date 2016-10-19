classdef GlobalFitPerformance < Figure
    % Figure comparing the performance of individual fits vs. a global fit.
    % A simulated dataset is fitted, both curve-by-curve individually, and
    % using a global fit. The global fit turns out to be much better at
    % extracting the correct values for the distance offset and force
    % offset (these are random offsets along the force and distance axes,
    % due to calibration errors and due to variations in microbead
    % diameter, respectively).

    properties

        % The distance offset added to the simulated data has a
        % Gaussian distribution with this sigma value (um).
        sigma_D      = 0.1;

        % The force offset added to the simulated data has a
        % Gaussian distribution with this sigma value (pN).
        sigma_F      = 0.5;

        % The simulated data is also stretched along the force
        % axis, to simulate the effect of errors in the force
        % calibration (calibration from PSD sensor voltage to
        % trap force, performed by measuring the power
        % spectral density of the Brownian motion of the bead
        % in the trap). Again, a Gaussian distribution with
        % this sigma value.
        sigma_deltaF = 0.01;

        % How many data curves to simulated.
        nCurves      = 100;

        % Handle to a function that generates a simulated dataset based
        % on the Odijk model.
        simulDataFun = @(Lc) simulatefddata(...
                                  'model',                 'odijk' ...
                                , 'modelParams',           [43 Lc 1766] ...
                                , 'dRange',                [12 17.5] ...
                                , 'fRange',                [0 30] ...
                                , 'noise',                 [0.25 0.005] ...
                                , 'nPoints',               500 ...
                                );

        % Contour length value (um) to use for the simulated dataset.
        Lc           = 16.5;

    end

    properties (SetAccess=private)

        % FdDataCollection holding the simulated F,d curves.
        simFd           = [];

        % Distance offset values used while simulating the
        % F,d curves (so we can compare the found values with
        % the true solutions).
        simulatedD0     = [];

        % Actual force offset values.
        simulatedF0     = [];

        % Actual force 'stretch' values (not used in the Figure).
        simulatedDeltaF = [];

        % Individual fit results.
        fitIndv         = [];

        % Global fit results.
        fitGlob         = [];

    end

    methods

        function [self] = GlobalFitPerformance(varargin)
            self = self@Figure(varargin{:});
        end

        function analyze(self)
            % ANALYZE Perform calculations.

            %% Generate simulated data.
            self.simulatedD0     = self.sigma_D      .* randn(self.nCurves,1);
            self.simulatedF0     = self.sigma_F      .* randn(self.nCurves,1);
            self.simulatedDeltaF = self.sigma_deltaF .* randn(self.nCurves,1) ...
                                   + 1.0;

            disp('Generating simulated data...');
            fdc = cell(self.nCurves,1);
            for i = 1:self.nCurves
                fd = self.simulDataFun(self.Lc);
                fd = fd.scale('f', self.simulatedDeltaF(i));
                fd = fd.shift('d', -self.simulatedD0(i));
                fd = fd.shift('f', self.simulatedF0(i));
                fdc{i} = fd;
            end
            self.simFd = FdDataCollection('skipduplicatecheck', fdc{:});


            %% Now fit it again.
            disp('Fitting simulated data...');
            self.fitIndv = fitfd(self.simFd, 'model', 'odijk-d0-f0', 'Lc', self.Lc);
            self.fitGlob = fitfdglobal(self.simFd, 'Lc', self.Lc);
        end

        function plot(self)
            % PLOT Create plot.

            plot@Figure(self);

            subplot(1,2,1);
            plotHandles = makePlot(self.simulatedD0, [self.fitIndv.d0], self.fitGlob.d0);
            title('d_0 (um)');

            legend(plotHandles, 'Individual fits', 'Global fit', 'Location', 'northwest');

            subplot(1,2,2);
            makePlot(self.simulatedF0, [self.fitIndv.F0], self.fitGlob.F0);
            title('F_0 (pN)');

            function [h] = makePlot(actual, foundIndv, foundGlob)
                lims = [min(actual) max(actual)];
                plot(lims, lims, '--', 'Color', [0.7 0.7 0.7]);
                hold('on');
                h1 = plot(actual, foundIndv, '.', 'Color', [0.8 0 0], 'MarkerSize', 10);
                h2 = plot(actual, foundGlob, '.', 'Color', [0 0.8 0], 'MarkerSize', 10);

                xlim(lims);
                ylim(lims);

                xlabel('Actual');
                ylabel('Found');

                h = [h1 h2];
            end

        end

    end

end
