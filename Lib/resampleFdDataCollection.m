function [resampledData] = resampleFdDataCollection(data)
% Returns an FdDataCollection with a resampled bootstrap dataset, created by
% resampling from the original data, with replacement.

nDatasets = data.length;
resampledData = FdDataCollection();
for k = datasample(1:nDatasets, nDatasets)
    resampledData.add('skipduplicatecheck', data.items{k});
end

end

