classdef MgVsFmax < Figure
    % The cut-off force "Fmax", vs. magnesium concentration.

    properties

        % TwlcVsMgAnalysis object. Can be left empty, in which case the
        % object will be constructed as necessary. If you already have
        % a TwlcVsMgAnalysis object containing analysis results, you can
        % also pass it in during construction of the figure:
        %
        %   >> f = Figure1('tvma', myTwlcVsMgAnalysis);
        tvma = [];

    end

    methods

        function [self] = MgVsFmax(varargin)
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

            Fmax = zeros(size(self.tvma.mgConcs));
            for i = 1:length(self.tvma.mgConcs)
                Fmax(i) = self.tvma.findCutOffForce(self.tvma.mgConcs(i));
            end

            plot(self.tvma.mgConcs, Fmax, '.b', 'MarkerSize', 20);

            xlabel('[Mg] (mM)');
            ylabel('F_{max} (pN)');
        end

    end

end
