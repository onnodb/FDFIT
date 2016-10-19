classdef DataAlignment < Figure
    % Comparing the alignment of F,d curves based on individual fitting,
    % vs. global fitting. The alignment based on the global fit performs
    % much better.

    properties

        % For this plot, only the first "nCurves" F,d curves (i.e., FdData
        % objects) from "data" are used. If we'd use all data, the
        % plot would get way too crowded.
        % (Note: "data" can be found in the parent class "Figure").
        nCurves               = 10;

        % For convenience, the first "nCurves" F,d curves that are
        % plotted are stored in this FdDataCollection, too.
        % You can also set this explicitly to an FdDataCollection yourself,
        % if you want to use your own data.
        dataUsed              = [];

        % The data from "dataUsed", aligned on the basis of individual
        % fits of each of the curves.
        dataAligned_IndivFits = [];

        % The data from "dataUsed", aligned on the basis of a global
        % fit of the curves.
        dataAligned_GlobFit   = [];

        % Needed for the second part of the alignment:
        % the distance region for the overstretching plateau.
        % Should be a vector [Dmin Dmax], with Dmin and Dmax both
        % in um.
        osPlateauRegion       = [18.5 20];

        % Matrix with xlim (1st row, um) and ylim (2nd row, pN) values.
        % Note: we're only interested in the eWLC regime for this
        % figure, so taking the upper limit of the Y axis as 30 pN.
        plotLims              = [12 17; -2 30];

    end

    % ------------------------------------------------------------------------

    methods

        function [self] = DataAlignment(varargin)
            self = self@Figure(varargin{:});
        end

        function analyze(self)
            % ANALYZE Perform calculations.

            % Fetch data to use: either pre-set data, or the first N
            % data curves in the dataset.
            if isempty(self.dataUsed)
                fprintf('Using first %d F,d curves in "data".\n', self.nCurves);
                self.dataUsed = FdDataCollection(self.data.items{1:self.nCurves});
            else
                disp('Using given FdDataCollection "dataUsed".');
            end

            % Initialize FdDataCollection objects.
            self.dataAligned_IndivFits = FdDataCollection();
            self.dataAligned_GlobFit   = FdDataCollection();


            %% Align the F,d curves using individual fits.
            disp('Performing alignment using individual fits...');
            for i = 1:self.nCurves
                f = fitfd(...
                            self.dataUsed.items{i} ...
                            , 'model',          'odijk-d0-f0' ...
                            , 'startParams',    [50 1500 0 0] ...
                            , 'Lc',             16.5 ...
                            );
                self.dataAligned_IndivFits.add(...
                        self.dataUsed.items{i}.shift('d', f.d0).shift('f', -f.F0) ...
                        );
            end


            %% Align the F,d curves, using global fit.
            disp('Performing alignment using global fit...');
            if isempty(self.dataAligned_GlobFit)
                disp('Aligning...');
                self.dataAligned_GlobFit = alignFdCurves(self.dataUsed, 'OSPlateauRegion', self.osPlateauRegion);
            else
                disp('Using pre-set aligned data.');
            end
        end

        function plot(self)
            % PLOT Create plot.

            plot@Figure(self);

            subplot(1,2,1);
            makePlots(self.dataAligned_IndivFits);
            title('(a) Alignment using individual fits', 'FontSize', 16);

            subplot(1,2,2);
            makePlots(self.dataAligned_GlobFit);
            title('(b) Alignment using global fit', 'FontSize', 16);

            function makePlots(fdc)
                for i = 1:fdc.length
                    fd = fdc.items{i}.subset('f', [-Inf 30]);
                    plot(fd.d, fd.f);
                    hold('on');
                end

                xlim(self.plotLims(1,:));
                ylim(self.plotLims(2,:));

                xlabel('Distance ({\mu}m)');
                ylabel('Force (pN)');
            end

        end

    end

end
