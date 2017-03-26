classdef BasicFdFitModel

    properties

        % Struct with fit parameter values.
        % Note: can only be initialized once.
        fitParams = struct();

        % Struct with fit parameter bounds. Each item is a length-2 vector
        % [min max].
        fitParamBounds = struct();

        % Struct with fixed parameter values.
        % Note: can only be initialized once.
        fixedParams = struct();

    end

    properties (SetAccess=protected)

        % Either 'd' or 'F'
        dependentVariable;

        % Either 'd' or 'F'
        independentVariable;

    end

    properties (Dependent)

        nFitParams;
        nFixedParams;

        hasFixedParams;

    end

    properties (Constant)

        default_Lp = 50;        % [nm]
        default_Lc = 16.5;      % [um]
        default_S  = 1500;      % [pN]
        default_kT = 4.11;      % [pN nm]
        default_Fc = 30.6;      % [pN]
        default_C  = 440;       % [pN nm^2]

    end

    % ------------------------------------------------------------------------

    methods

        function [val] = get.hasFixedParams(self)
            val = ~isempty(fieldnames(self.fixedParams));
        end

        function [val] = get.nFitParams(self)
            val = length(fieldnames(self.fitParams));
        end

        function [val] = get.nFixedParams(self)
            val = length(fieldnames(self.fixedParams));
        end

        function [self] = set.fitParams(self, val)
            if isempty(fieldnames(self.fitParams)) && isempty(fieldnames(val))
                % OK
            elseif isempty(fieldnames(self.fitParams))
                % Initialization is OK
                self.fitParams = val;
                pb = struct();
                fn = fieldnames(val);
                for i = 1:length(fn)
                    pb.(fn{i}) = [-Inf +Inf];
                end
                self.fitParamBounds = pb;
            elseif isstruct(val)
                self.validateStructAssignment(val, self.fitParams);
                self.fitParams = val;
            elseif isnumeric(val)
                self.validateNumAssignment(val, self.fitParams);
                fn = fieldnames(self.fitParams);
                for i = 1:length(fn)
                    self.fitParams.(fn{i}) = val(i);
                end
            else
                error('Invalid value for "fitParams".');
            end
        end

        function [self] = set.fitParamBounds(self, val)
            if isempty(fieldnames(self.fitParamBounds)) && isempty(fieldnames(val))
                % OK
            elseif isempty(fieldnames(self.fitParamBounds))
                % Initialization is OK
                self.fitParamBounds = val;
            elseif isstruct(val)
                self.validateStructAssignment(val, self.fitParamBounds);
                fn = fieldnames(val);
                for i = 1:length(fn)
                    v = val.(fn{i});
                    if ~isnumeric(v) || ~isvector(v) || length(v) ~= 2
                        error('Invalid bounds vector for "fitParamBounds.%s".', fn{i});
                    end
                end
                self.fitParamBounds = val;
            elseif isnumeric(val)
                if ~ismatrix(val) || size(val,1) ~= length(fieldnames(self.fitParamBounds)) ...
                        || size(val,2) ~= 2
                    error('Invalid numeric assignment for "fitParamBounds".');
                end
                fn = fieldnames(self.fitParamBounds);
                for i = 1:length(fn)
                    self.fitParamBounds.(fn{i}) = val(i,:);
                end
            else
                error('Invalid value for "fitParamBounds".');
            end
        end

        function [self] = set.fixedParams(self, val)
            if isempty(fieldnames(self.fixedParams)) && isempty(fieldnames(val))
                % OK
            elseif isempty(fieldnames(self.fixedParams))
                % Initialization is OK
                self.fixedParams = val;
            elseif isstruct(val)
                self.validateStructAssignment(val, self.fixedParams);
                self.fixedParams = val;
            elseif isnumeric(val)
                self.validateNumAssignment(val, self.fixedParams);
                fn = fieldnames(self.fixedParams);
                for i = 1:length(fn)
                    self.fixedParams.(fn{i}) = val(i);
                end
            else
                error('Invalid value for "fixedParams".');
            end
        end

    end

    methods (Access=private)

        function validateStructAssignment(~, val, oldval)
            if ~isempty(setdiff(fieldnames(val), fieldnames(oldval)))
                error('Invalid struct assignment.');
            end
        end

        function validateNumAssignment(~, val, oldval)
            if ~isvector(val) || length(val) ~= length(fieldnames(oldval))
                error('Invalid numeric assignment.');
            end
        end

    end

end
