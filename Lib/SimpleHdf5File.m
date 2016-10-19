classdef SimpleHdf5File < handle
    % SIMPLEHDF5FILE Simple wrapper around MATLAB's low-level HDF5 functions.

    properties (SetAccess = protected)
        fileId      = [];
        filename    = [];
        writeAccess = false;
    end

    methods

        function [self] = SimpleHdf5File(filename, accessMode)
            % Open or create an HDF5 file.
            %
            % INPUT:
            % filename = name of the hdf5 file. An '.h5' extension is
            %       added automatically if 'filename' doesn't include any
            %       extension.
            %       NOTE: On UNIX, the tilde (~) is NOT expanded.
            % accessMode = r|w|a|r+|w+|a+ (default: 'r'; as in 'fopen')
            %       'r':  Open file for reading.
            %       'w':  Open or create new file for writing. Discard existing
            %             contents, if any.
            %       'a':  Open or create new file for writing. If file exists,
            %             open for reading/writing ("append" mode).
            %       'r+': Open file for reading and writing.
            %       'w+': (effectively the same as 'w')
            %       'a+': (effectively the same as 'a')
            %
            % SEE ALSO:
            % fopen

            if nargin < 2
                accessMode = 'r';
            end

            [~, ~, fileExt] = fileparts(filename);
            if isempty(fileExt)
                filename = [filename '.h5'];
            end

            switch accessMode
                case 'r'
                    self.fileId = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
                case {'w','w+'}
                    self.fileId = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
                    self.writeAccess = true;
                case {'a','a+'}
                    if exist(filename, 'file')
                        self.fileId = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
                    else
                        self.fileId = H5F.create(filename, 'H5F_ACC_EXCL', 'H5P_DEFAULT', 'H5P_DEFAULT');
                    end
                    self.writeAccess = true;
                case 'r+'
                    self.fileId = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
                    self.writeAccess = true;
                otherwise
                    error('SIMPLEHDF5FILE:invalidArgument', 'Invalid access mode "%s".', accessMode);
            end

            self.filename   = filename;
        end

        function delete(self)
            H5F.close(self.fileId);
            self.fileId = [];
        end

        function [attr] = getAttributes(self, path)
            % GETATTRIBUTES Get the attributes attached to an object.
            %
            % INPUT:
            % path = path to the object to read the attributes from.
            %       Use '/' to get attributes attached to the root.
            %
            % OUTPUT:
            % attr = a scalar struct, with field names indicating attribute
            %       names, and field values giving the corresponding values.

            attr = struct();

            objId = H5O.open(self.fileId, path, 'H5P_DEFAULT');
            objId_cleanup = onCleanup(@() H5O.close(objId));
            info = H5O.get_info(objId);

            for i = 0:info.num_attrs-1
                attrId = H5A.open_by_idx(self.fileId, path, 'H5_INDEX_CRT_ORDER', ...
                                         'H5_ITER_INC', i, 'H5P_DEFAULT', 'H5P_DEFAULT');
                attrName = H5A.get_name(attrId);
                attrValue = self.readAttribute(attrId);
                H5A.close(attrId);

                attr.(attrName) = attrValue;
            end
        end

        function [data] = readData(self, datasetPath)
            % READDATA Read data from a dataset.
            %
            % INPUT:
            % datasetPath = string representing the (full) path of the dataset
            %       to write. Any intermediate (nested) groups are created
            %       automatically.
            %
            % OUTPUT:
            % data = any MATLAB value, e.g., a matrix, string, etc.
            %       If the datatype is not supported, an exception will be
            %       raised.

            dsetId = H5D.open(self.fileId, datasetPath);
            dsetId_cleanup = onCleanup(@() H5D.close(dsetId));

            spaceId = H5D.get_space(dsetId);
            spaceId_cleanup = onCleanup(@() H5S.close(spaceId));

            dataTypeId = H5D.get_type(dsetId);
            dataTypeId_cleanup = onCleanup(@() H5T.close(dataTypeId));

            if H5T.get_class(dataTypeId) == H5ML.get_constant_value('H5T_COMPOUND')
                data = self.readCompoundData(dsetId, dataTypeId, spaceId);
            else
                data = H5D.read(dsetId);
            end
            data = self.postProcessData(data, dataTypeId);
        end

        function setAttributes(self, path, attr)
            % SETATTRIBUTES Attach a set of attributes to an object.
            %
            % INPUT:
            % path = path to the object to attach the attributes to.
            %       Use '/' to attach attributes to the root.
            % attr = a scalar struct, with field names indicating attribute
            %       names, and field values giving the corresponding values.

            objId = H5O.open(self.fileId, path, 'H5P_DEFAULT');
            objId_cleanup = onCleanup(@() H5O.close(objId));

            fn = fieldnames(attr);
            for i = 1:length(fn)
                self.writeAttribute(objId, fn{i}, attr.(fn{i}));
            end
        end

        function writeData(self, datasetPath, data)
            % WRITEDATA Write data to a dataset.
            %
            % INPUT:
            % datasetPath = string representing the (full) path of the dataset
            %       to write. Any intermediate (nested) groups are created
            %       automatically.
            % data = any MATLAB value, e.g., a matrix, string, etc.
            %       If the datatype is not supported, an exception will be
            %       raised.

            %% Validation.
            if ~self.writeAccess
                error('SIMPLEHDF5FILE:readOnlyMode', 'Cannot write data: file open in read-only mode.');
            end

            %% Prep data if needed.
            if isstruct(data)
                data = self.convertStringsInStructToCellstrings(data);
            end

            %% Do the magic.
            if self.datasetExists(datasetPath)
                % Rewrite existing dataset; but this only works if the
                % dimensions match up...
                self.writeDataToExistingDataset(datasetPath, data);
            else
                % Create new dataset, and write the data to it.
                self.writeDataToNewDataset(datasetPath, data);
            end

        end

        function [b] = datasetExists(self, datasetPath)
            try
                dsId = H5D.open(self.fileId, datasetPath);
                H5D.close(dsId);
                b = true;
            catch
                b = false;
            end
        end

    end

    methods (Access = private)

        function [b] = attributeExists(self, objId, attrName)
            try
                attrId = H5A.open(objId, attrName, 'H5P_DEFAULT');
                H5A.close(attrId);
                b = true;
            catch
                b = false;
            end
        end

        function [data] = postProcessData(self, data, dataTypeId)
            switch H5T.get_class(dataTypeId)
                case H5ML.get_constant_value('H5T_STRING')
                    if isvector(data)
                        % MATLAB strings are row vectors, now column vectors.
                        data = data';
                    end
                case H5ML.get_constant_value('H5T_COMPOUND')
                    data = self.convertCellstringsInStructToStrings(data);
            end
        end

        function [v] = readAttribute(self, attrId)
            v = H5A.read(attrId, 'H5ML_DEFAULT');

            spaceId = H5A.get_space(attrId);
            spaceId_cleanup = onCleanup(@() H5S.close(spaceId));

            dataTypeId = H5A.get_type(attrId);
            dataTypeId_cleanup = onCleanup(@() H5T.close(dataTypeId));

            v = self.postProcessData(v, dataTypeId);
        end

        function writeAttribute(self, objId, attrName, attrValue)
            %% Get datatype.
            dataTypeId = self.dataToHdfDataType(attrValue);
            dataTypeId_cleanup = onCleanup(@() H5T.close(dataTypeId));

            %% Create dataspace (using datatype + dims).
            spaceId = self.dataToDataSpace(attrValue);
            spaceId_cleanup = onCleanup(@() H5S.close(spaceId));

            %% Create, or if necessary, re-create the attribute.
            acpl = H5P.create('H5P_ATTRIBUTE_CREATE');
            acpl_cleanup = onCleanup(@() H5P.close(acpl));

            if self.attributeExists(objId, attrName)
                % Make sure the datatype & dataspace match up, otherwise
                % we'll have to re-create the attribute.
                attrId = H5A.open(objId, attrName, 'H5P_DEFAULT');

                curDataType  = H5A.get_type(attrId);
                curDataType_cleanup = onCleanup(@() H5T.close(curDataType));
                curSpace     = H5A.get_space(attrId);
                curSpace_cleanup = onCleanup(@() H5S.close(curSpace));
                [~, curDims] = H5S.get_simple_extent_dims(curSpace);

                if (~H5T.equal(dataTypeId, curDataType)) || (~isequal(fliplr(curDims), size(attrValue)))
                    % We have to delete the attribute & recreate it.
                    H5A.close(attrId);
                    H5A.delete(objId, attrName);

                    attrId = H5A.create(objId, attrName, dataTypeId, spaceId, acpl);
                end
            else
                attrId = H5A.create(objId, attrName, dataTypeId, spaceId, acpl);
            end

            H5A.write(attrId, dataTypeId, attrValue);
            H5A.close(attrId);
        end

        function writeDataToExistingDataset(self, datasetPath, data)
            error('SIMPLEHDF5FILE:operationNotSupported', ...
                  'Writing to an existing dataset is not currently supported.');
        end

        function writeDataToNewDataset(self, datasetPath, data)
            %% Get datatype.
            dataTypeId = self.dataToHdfDataType(data);
            dataTypeId_cleanup = onCleanup(@() H5T.close(dataTypeId));

            %% Create dataspace (using datatype + dims).
            spaceId = self.dataToDataSpace(data);
            spaceId_cleanup = onCleanup(@() H5S.close(spaceId));

            %% Create dataset.
            % Indicate we want to automatically create any intermediate
            % (nested) groups whenever necessary.
            lcpl = H5P.create('H5P_LINK_CREATE');
            lcpl_cleanup = onCleanup(@() H5P.close(lcpl));
            H5P.set_create_intermediate_group(lcpl, 1);

            % Create dataset.
            dsetId = H5D.create(self.fileId, datasetPath, dataTypeId, spaceId, ...
                                lcpl, 'H5P_DEFAULT', 'H5P_DEFAULT');
            dsetId_cleanup = onCleanup(@() H5D.close(dsetId));

            %% Write data.
            if isstruct(data)
                memSpaceId = H5S.create_simple(1, 1, []);
                for i = 1:length(data)
                    H5S.select_hyperslab(spaceId, 'H5S_SELECT_SET', i-1, [], [], []);
                    H5D.write(dsetId, 'H5ML_DEFAULT', memSpaceId, spaceId, 'H5P_DEFAULT', data(i));
                end
                H5S.close(memSpaceId);
            else
                H5D.write(dsetId, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data);
            end
        end

        function [spaceId] = dataToDataSpace(~, data)
            if isempty(data)
                spaceId = H5S.create('H5S_NULL');
            elseif ischar(data)
                spaceId = H5S.create('H5S_SCALAR');
            elseif iscell(data) && ~isempty(data) && ischar(data{1}) && isvector(data)
                spaceId = H5S.create_simple(1, length(data), length(data));
            elseif isstruct(data) && isvector(data)
                spaceId = H5S.create_simple(1, length(data), length(data));
            else
                spaceId = H5S.create_simple(ndims(data), fliplr(size(data)), fliplr(size(data)));
            end
        end

        function [dataTypeId] = dataToHdfDataType(self, data)
            switch class(data)
                case 'single'
                    dataTypeId = H5T.copy('H5T_NATIVE_FLOAT');
                case 'double'
                    dataTypeId = H5T.copy('H5T_NATIVE_DOUBLE');
                case 'int8'
                    dataTypeId = H5T.copy('H5T_NATIVE_CHAR');
                case 'int16'
                    dataTypeId = H5T.copy('H5T_NATIVE_SHORT');
                case 'int32'
                    dataTypeId = H5T.copy('H5T_NATIVE_INT');
                case 'int64'
                    dataTypeId = H5T.copy('H5T_NATIVE_INT64');
                case 'uint8'
                    dataTypeId = H5T.copy('H5T_NATIVE_UCHAR');
                case 'uint16'
                    dataTypeId = H5T.copy('H5T_NATIVE_USHORT');
                case 'uint32'
                    dataTypeId = H5T.copy('H5T_NATIVE_UINT');
                case 'uint64'
                    dataTypeId = H5T.copy('H5T_NATIVE_UINT64');
                case 'logical'
                    dataTypeId = H5T.copy('H5T_NATIVE_CHAR');
                case 'char'
                    if ~isempty(data) && ~isrow(data)
                        error('SIMPLEHDF5FILE:invalidDataType', ...
                              'Only single strings (row vectors of type "char") are supported.');
                    end

                    dataTypeId = H5T.copy('H5T_C_S1');
                    if ~isempty(data)
                        H5T.set_size(dataTypeId, numel(data));
                    end
                    H5T.set_strpad(dataTypeId, 'H5T_STR_NULLTERM');
                case 'cell'
                    if isempty(data)
                        error('SIMPLEHDF5FILE:unsupportedDataType', ...
                              'Empty cell arrays are not supported.');
                    end
                    if isvector(data) && ischar(data{1})
                        % We do support cell arrays of strings.
                        dataTypeId = H5T.vlen_create(H5T.copy('H5T_C_S1'));
                    else
                        error('SIMPLEHDF5FILE:unsupportedDataType', ...
                              'This type of cell array is not supported.');
                    end
                case 'struct'
                    if ~isscalar(data) && ~iscolumn(data)
                        error('SIMPLEHDF5FILE:invalidDataType', ...
                            'Only struct scalars and column vectors are supported.');
                    end

                    dataTypeId = self.createCompoundDataType(data);
                otherwise
                    error('SIMPLEHDF5FILE:invalidDataType', ...
                          'Invalid datatype "%s".', class(data));
            end
        end

        function [dataTypeId] = createCompoundDataType(self, data)
            assert(isstruct(data));

            fieldNames = fieldnames(data);
            nFields    = length(fieldNames);

            %% Collect data types and corresponding data sizes.
            compDataTypes  = cell(nFields,1);
            compFieldSizes = zeros(nFields,1);

            for i = 1:nFields
                curFieldValue = data(1).(fieldNames{i});

                if isnumeric(curFieldValue) && ~isscalar(curFieldValue)
                    error('SIMPLEHDF5FILE:invalidDataType', ...
                          'Only structs with scalar fields are supported.');
                end
                compDataTypes{i} = self.dataToHdfDataType(curFieldValue);

                compFieldSizes(i) = H5T.get_size(compDataTypes{i});
            end

            %% Calculate field offsets (in bytes).
            compFieldOffsets = zeros(nFields,1);
            compFieldOffsets(2:end) = cumsum(compFieldSizes(1:end-1));

            %% Create compound data type.
            dataTypeId = H5T.create('H5T_COMPOUND', sum(compFieldSizes));
            for i = 1:nFields
                H5T.insert(dataTypeId, fieldNames{i}, compFieldOffsets(i), compDataTypes{i});
            end
        end

        function [data] = readCompoundData(self, dsetId, dataTypeId, spaceId)
            if H5S.get_simple_extent_ndims(spaceId) ~= 1
                error('SIMPLEHDF5FILE:cannotReadHighDimCompoundData', ...
                      'Reading in compound data of rank > 1 is not currently supported.');
            end

            nFields            = H5T.get_nmembers(dataTypeId);
            compDataTypes      = cell(nFields,1);
            compDataFieldNames = cell(nFields,1);

            for i = 1:nFields
                compDataTypes{i}      = H5T.get_member_type(dataTypeId, i-1);
                compDataFieldNames{i} = H5T.get_member_name(dataTypeId, i-1);
            end

            nPoints = H5S.get_simple_extent_npoints(spaceId);

            memSpaceId = H5S.create_simple(1, 1, []);
            for i = 1:nPoints
                H5S.select_hyperslab(spaceId, 'H5S_SELECT_SET', i-1, [], [], []);
                data(i) = H5D.read(dsetId, 'H5ML_DEFAULT', memSpaceId, spaceId, 'H5P_DEFAULT');
                for j = 1:nFields
                    data(i).(compDataFieldNames{j}) = self.postProcessData(data(i).(compDataFieldNames{j}), compDataTypes{j});
                end
            end
            H5S.close(memSpaceId);

            data = data';
        end

        function [data] = convertStringsInStructToCellstrings(self, data)
            assert(isstruct(data));

            fn = fieldnames(data);
            for i = 1:length(data)
                for j = 1:length(fn)
                    if ischar(data(i).(fn{j}))
                        data(i).(fn{j}) = {data(i).(fn{j})};    % make into cellstring
                    end
                end
            end
        end

        function [data] = convertCellstringsInStructToStrings(self, data)
            assert(isstruct(data));

            fn = fieldnames(data);
            for i = 1:length(data)
                for j = 1:length(fn)
                    if iscellstr(data(i).(fn{j})) && isscalar(data(i).(fn{j}))
                        data(i).(fn{j}) = data(i).(fn{j}){1};    % back into character array
                    end
                end
            end
        end

    end

end
