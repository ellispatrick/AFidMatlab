function [maskAF, im1AFRemoved, im2AFRemoved, kBest] = afIdentifier(im1, im2, bw, varargin)
% Inputs:
% im1 - image 1
% im2 - image 2
% bw  - mask of objects in the image to measure
%
% Name-value pair arguments:
% 'Corr'  - correlation cut-off for autofluorescence; all objects above
%         - this cut-off within the identified autofluorescence cluster is
%         - classified as autofluorescence
%         - Values taken: values between -1 and 1
%
% 'kAuto' - estimate the number of clusters between 3 and a user provided
%         - number of clusters defined by 'k'
%         - Values taken: 0 = do not estimate number of clusters
%                         1 = estimate number of clusters
%
% 'k'     - if kAuto = 0: number of clusters used when performing k-means
%         - if kAuto = 1: maximum number of clusters for estimating an
%                         optimal k when performing k-means
%         - Values taken: positive integers; if k = 1, no k-means is
%                         performed
%
% Outputs:
% maskAF       - mask of autofluorescent objects
% im1AFRemoved - im1 with autofluorescent objects set to 0
% im2AFRemoved - im2 with autofluorescent objects set to 0
% kBest        - estimated value of the optimum k if automated (the value
%              - of k used if not automated)

%% Parameters
p = inputParser;

defaultK = 1;
defaultCorr = -1;
defaultKAuto = 0;

addRequired(p, 'im1', @ismatrix);
addRequired(p, 'im2', @ismatrix);
addRequired(p, 'bw', @ismatrix);

addParameter(p, 'k', defaultK, @isnumeric);
addParameter(p, 'Corr', defaultCorr, @isnumeric);
addParameter(p, 'kAuto', defaultKAuto, @isnumeric);

parse(p, im1, im2, bw, varargin{:});

% Read in images
im1 = p.Results.im1;
im2 = p.Results.im2;
bw = p.Results.bw;

if ~isequal(size(im1), size(im2)) && ~isequal(size(im1), size(bw))
    error('Image sizes are not equal.');
end

% Read in parameters
kUser = p.Results.k;

if kUser < 1 || floor(kUser) ~= kUser
    error('k must be a positive integer.');
end

corrCutoff = p.Results.Corr;

if corrCutoff < -1 || corrCutoff > 1
    error('Corr must be between -1 and 1.');
end

kAuto = p.Results.kAuto;

if kAuto ~= 0 && kAuto ~= 1 || floor(kAuto) ~= kAuto
    error('k must be a positive integer.');
elseif kAuto == 1 && kUser <= 3
    error('k for automatically estimating k must be greater than 3.');
end

%% Classify autofluorescence
im1PixelsStruct = regionprops(bw, im1, 'PixelValues');
objCount = length(im1PixelsStruct);

im2PixelsStruct = regionprops(bw, im2, 'PixelValues');

pixelsStruct = regionprops(bw, 'PixelIdxList');

% Correlation
corrVals = zeros(1, objCount);

for i = 1:objCount
    im1PixelVals = double(im1PixelsStruct(i).PixelValues);
    im2PixelVals = double(im2PixelsStruct(i).PixelValues);
    
    corrVals(i) = corr(im1PixelVals, im2PixelVals);
end

corrAF = corrVals > corrCutoff;

% K-means
if kUser > 1
    
    % Additional texture measures
    std1Vals = zeros(1, objCount);
    std2Vals = zeros(1, objCount);
    kurt1Vals = zeros(1, objCount);
    kurt2Vals = zeros(1, objCount);
    
    for i = 1:objCount
        im1PixelVals = double(im1PixelsStruct(i).PixelValues);
        im2PixelVals = double(im2PixelsStruct(i).PixelValues);
        
        std1Vals(i) = std(im1PixelVals);
        std2Vals(i) = std(im2PixelVals);
        kurt1Vals(i) = kurtosis(im1PixelVals);
        kurt2Vals(i) = kurtosis(im2PixelVals);
    end
    
    % Transform texture measures for clustering
    corrValsTform = atanh(corrVals);
    corrValsTform = corrValsTform/std(corrValsTform);
    std1ValsTform = log(std1Vals);
    std1ValsTform = std1ValsTform/std(std1ValsTform);
    std2ValsTform = log(std2Vals);
    std2ValsTform = std2ValsTform/std(std2ValsTform);
    kurt1ValsTform = log(kurt1Vals);
    kurt1ValsTform = kurt1ValsTform/std(kurt1ValsTform);
    kurt2ValsTform = log(kurt2Vals);
    kurt2ValsTform = kurt2ValsTform/std(kurt2ValsTform);
    
    toCluster = [corrValsTform' std1ValsTform' std2ValsTform' ...
        kurt1ValsTform' kurt2ValsTform'];
    
    if kAuto % Automated
        kMax = kUser;
        
        % Initialise arrays
        idx = zeros(kMax, objCount);
        
        corrClust = zeros(1, kMax);
        
        statVals = zeros(1,kMax);
        
        afSize = zeros(1,kMax);
        
        % Perform k-means
        for k = 1:kMax
            idx(k,:) = transpose(kmeans(toCluster, k, 'MaxIter', 10000));
            
            % Identify top two clusters
            corrMean = splitapply(@mean, corrValsTform, idx(k,:));
            
            [~, corrClust(k)] = max(corrMean);
            corrMean(corrClust(k)) = -Inf;
            [~, secondClust] = max(corrMean);
            
            afSize(k) = sum(idx(k,:) == corrClust(k));
            
            % Measure t test values
            corrVals1 = corrValsTform(idx(k,:) == corrClust(k));
            corrVals2 = corrValsTform(idx(k,:) == secondClust);
            
            [~, ~, ~, stats] = ttest2(corrVals1, corrVals2);
            
            statVals(k) = stats.tstat;
        end
        
        %%% Perform elbow test
        x1 = 3;
        x2 = kMax;
        y1 = statVals(3);
        y2 = statVals(kMax);
        
        [~, kBest] = max((y2-y1)*(x1:x2) - (x2-x1)*statVals(x1:x2) + x2*y1-y2*x1);
        
        kBest = kBest + 2;
        
        clustAF = (idx(kBest,:) == corrClust(kBest));
        
    else % Non-automated
        idx = kmeans(toCluster, kUser, 'MaxIter', 10000);
        corrMean = splitapply(@mean, corrValsTform, idx);
        [~, afIdx] = max(corrMean);
        clustAF = (idx == afIdx);
        
        
        kBest = kUser;
    end
else
    clustAF = ones(1, objCount);
    
end

afIdx = corrAF & clustAF;

%% Remove autofluorescence
maskAF = bw;
 
for i = 1:objCount
    pixels = pixelsStruct(i).PixelIdxList;
    
    if ~afIdx(i)
        maskAF(pixels) = 0;
    end
end
 
im1AFRemoved = im1;
im1AFRemoved(maskAF==0) = 0;
 
im2AFRemoved = im2;
im2AFRemoved(maskAF==0) = 0;

