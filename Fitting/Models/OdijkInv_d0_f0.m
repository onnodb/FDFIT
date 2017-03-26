classdef OdijkInv_d0_f0 < BasicFdFitModel

    methods

        function [self] = OdijkInv_d0_f0()
            self.dependentVariable   = 'F';
            self.independentVariable = 'd';

            self.fitParams = struct(...
                                  'Lp',     BasicFdFitModel.default_Lp ...
                                , 'S',      BasicFdFitModel.default_S ...
                                , 'd0',     0 ...
                                , 'F0',     0 ...
                                );

            self.fixedParams = struct(...
                                  'Lc',     BasicFdFitModel.default_Lc ...
                                );

            self.fitParamBounds.F0 = [-20 20];
        end

        function [ftype] = getFitType(self)
            ftype = fittype(...
                        @(Lp, S, d0, F0, x) fOdijkInv_d0_f0(...
                                                x, Lp, self.fixedParams.Lc, ...
                                                S, d0, F0) ...
                        );
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
