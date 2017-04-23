classdef Odijk_Lc < BasicFdFitModel

    methods

        function [self] = Odijk_Lc()
            self.dependentVariable   = 'd';
            self.independentVariable = 'F';

            self.fitParams = struct(...
                                  'Lp',     BasicFdFitModel.default_Lp ...
                                , 'S',      BasicFdFitModel.default_S ...
                                , 'Lc',     BasicFdFitModel.default_Lc ...
                                );
        end

        function [fun] = getFitFun(self)
            fun = ...
                @(Lp, S, Lc, x) fOdijk_d0_f0(x, Lp, Lc, S, 0, 0) ...
                ;
        end

        function [fd_fit] = validateData(self, fd)
            fd_fit = fd;
            if any(fd.f > 30)
                warning('Only using data up to 30 pN.');
                fd_fit = fd.subset('f', [-Inf 30]);
            end
        end

    end

end
