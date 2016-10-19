classdef Figure < handle
    % FIGURE Base class encapsulating a figure for in the paper / SI.

    properties

        % FdDataCollection with all the data for the paper.
        data                = [];

        % Handle to the figure window.
        fig                 = [];

        axes_watermark      = [];

        % Set to true to automatically export a PNG of the figure.
        autoExport          = false;

    end % properties

    % ------------------------------------------------------------------------

    methods

        function [self] = Figure(varargin)
            % FIGURE Constructor.
            %
            % SYNTAX:
            % f = Figure('key', value, ...)
            %
            % KEY-VALUE PAIR ARGUMENTS:
            % All class properties can be initialized using the key-value pair
            % system. (Note: these are case-sensitive).
            %
            % EXAMPLE:
            % >> myfdc = FdDataCollection();
            % >> f = Figure('data', myfdc, 'autoExport');
            %
            % NOTE: If no value for "data" is given in the key-value pairs,
            % an attempt is made to automatically load the data from
            % a variable called "data" in the global workspace.
            %
            % NOTE: Automatically calls the implementation class's "analyze"
            % and "plot" methods, and, if the "autoExport" flag is set,
            % "export".

            % Parse key-value pair arguments
            if ~isempty(varargin)
                if isa(varargin{1}, 'FdDataCollection')
                    self.data = FdDataCollection(varargin{1});
                    varargin = varargin(2:end);
                end

                parseClassArgs(varargin, self);
            end

            % Get data from base workspace, if needed
            if isempty(self.data)
                if evalin('base', 'exist(''data'', ''var'') ~= 1')
                    error('No data given, and also failed to find data in global workspace');
                end
                self.data = FdDataCollection(evalin('base', 'data'));
            end

            % Auto-run
            try
                self.analyze();
                self.plot();
                if self.autoExport
                    self.export();
                end
            catch errorObj
                disp(getReport(errorObj, 'extended'));
                return
            end
        end

        function export(self, filename)
            % EXPORT Export the figure to a PNG file.
            %
            % INPUT:
            % filename = filename of the PNG file to export to.
            %       (Default: a file "Figure.png" in the system's temp
            %       directory).

            if nargin < 2
                filename = fullfile(tempdir(), [class(self) '.png']);
            end

            print(self.fig, filename, '-painters', '-dpng', '-r300');
        end

        function plot(self)
            % PLOT Creates the actual plot in the Figure.
            %
            % Subclasses should implement this function to create the actual
            % plot.

            self.fig = figure('Name', class(self), 'Color', [1 1 1]);
        end

    end % methods

    % ------------------------------------------------------------------------

    methods (Abstract)

        analyze(self);
            % ANALYZE Perform any calculations that are needed for the plot.
            %
            % Subclasses should implement this function to create the actual
            % plot.

    end % methods (Abstract)

    % ------------------------------------------------------------------------

    methods (Access = protected)

        function setupAxes(self)
            % Helper function that configures the current axes; sets
            % things like font size, etc.

            set(gca(), 'FontSize', 14);
        end

    end % methods (protected)

end
