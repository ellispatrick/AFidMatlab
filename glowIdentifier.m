function [glow1, glow2, im1GlowRemoved, im2GlowRemoved] = glowIdentifier(im1, im2, bw, varargin)
% Inputs:
% im1 - image 1 (before autofluorescence removal)
% im2 - image 2 (before autofluorescence removal)
% bw  - mask of autofluorescent objects identified by afRemove.m
%
% Name-value pair arguments:
% 'TraceSensitivity'  - number of steps taken by tracing algorithm before
%                     - dropping a point of expansion
%                     - Values taken: positive integers greater than 1
%
% 'Sigma'             - sigma for Gaussian blurring the image for measuring
%                     - pixel values
%                     - Values taken: values greater than or equal to 0 
%                                     (if Sigma = 0, images are not blurred)
%
% Outputs:
% glow1          - mask of glow in im1
% glow2          - mask of glow in im2
% im1GlowRemoved - im1 with glow set to 0
% im2GlowRemoved - im2 with glow set to 0
%% Parameters
p = inputParser;

defaultTraceSensitivity = 20;
defaultSigma = 0;

addRequired(p, 'im1', @ismatrix);
addRequired(p, 'im2', @ismatrix);
addRequired(p, 'bw', @ismatrix);

addParameter(p, 'TraceSensitivity', defaultTraceSensitivity, @isnumeric);
addParameter(p, 'Sigma', defaultSigma, @isnumeric);

parse(p, im1, im2, bw, varargin{:});

% Read in images
im1 = p.Results.im1;
im2 = p.Results.im2;
bw = p.Results.bw;

if ~isequal(size(im1), size(im2)) && ~isequal(size(im1), size(bw))
    error('Image sizes are not equal.');
end

% Read in parameters
traceSensitivity = p.Results.TraceSensitivity;

if traceSensitivity <= 1 || floor(traceSensitivity) ~= traceSensitivity
    error('Trace sensitivity should be an integer greater than 1.');
end

sigma = p.Results.Sigma;

if sigma < 0
    error('Sigma should be greater than or 0.');
end

%% Identify points of expansion
skelAF = bwmorph(bw, 'skel', Inf);

% Loop through entire image to identify end nodes (pixels in the skeleton
% which only have 1 neighbouring pixel); these end nodes serve as seeds for
% a tracing algorithm

[m, n] = size(im1);

endNodes = false(m,n);

for r = 1:m
    for c = 1:n
        if skelAF(r,c) == 1
            neighbours = zeros(1,8);
            
            rNeighbour = [r-1 r r+1 r-1 r+1 r-1 r r+1];
            cNeighbour = [c-1 c-1 c-1 c c c+1 c+1 c+1];
            
            for i = 1:length(rNeighbour)
                try
                    neighbours(i) = skelAF(rNeighbour(i),cNeighbour(i));
                catch
                    neighbours(i) = 0;
                end
            end
            
            if sum(neighbours) < 2
                endNodes(r,c) = 1;
            end
        end
    end
end

% Annular objects might not have end nodes - this ensures that a seed is
% defined for these objects
withNode = cell2mat(struct2cell(regionprops(skelAF, endNodes, 'MaxIntensity')));

skeletonPixelIdx = regionprops(skelAF, 'PixelIdxList');

for i = 1:length(skeletonPixelIdx)
    if ~withNode(i)
        endNodes(skeletonPixelIdx(i).PixelIdxList(1)) = 1;
    end
end

% Tracing algorithm
trace = endNodes;

expPoints = endNodes;

done = false(m,n);

traceCount = 0;

while 1
    traceCount = traceCount+1;
    
    newTrace = false(m,n);
    
    for r = 1:m
        for c = 1:n
            if trace(r,c)
                
                done(r,c) = 1;
                
                rNeighbour = [r-1 r r+1 r-1 r+1 r-1 r r+1];
                cNeighbour = [c-1 c-1 c-1 c c c+1 c+1 c+1];
                
                for i = 1:length(rNeighbour)
                    r2 = rNeighbour(i);
                    c2 = cNeighbour(i);
                    try
                        if (r2 > 0 && r2 <= m) && (c2 > 0 && c2 <= n) ...
                                && skelAF(r2,c2) ...
                                && ~done(r2,c2) ...
                                && ~trace(r2,c2) ...
                                && ~newTrace(r2,c2)
                            newTrace(r2,c2) = 1;
                            
                            if ~mod(traceCount,traceSensitivity)
                                expPoints(r2,c2) = 1;
                            end
                        end
                    catch
                        newTrace(r2,c2) = 0;
                    end
                end
                
            end
        end
    end
    
    trace = newTrace;
    
    if ~sum(trace(:))
        break;
    end
end

%%% Gaussian blur images
if sigma > 0
    im1Blurred = imgaussfilt(im1,sigma);
    im2Blurred = imgaussfilt(im2,sigma);
else
    im1Blurred = im1;
    im2Blurred = im2;
end

%% Identify body
extendDist = 60;
extendWidth = 1;

th = acos(1-(extendWidth/extendDist)^2/2);

minSteps = 3;
maxSteps = 30;

glow1 = false(m,n);
glow2 = false(m,n);

for r = 1:m
    for c = 1:n
        if expPoints(r,c)
            thMultiple = 1;
            thMeasure = th;
            
            while thMeasure < 2*pi
                thMeasure = th * thMultiple;
                xStepSize = cos(thMeasure);
                yStepSize = sin(thMeasure);
                
                extendLength = 1;
                extendSteps = 0;
                
                done1 = 0;
                done2 = 0;
                
                while 1
                    xMeasure = floor(r + extendLength * xStepSize);
                    yMeasure = floor(c + extendLength * yStepSize);
                    
                    xCompare = floor(r + (extendLength-1) * xStepSize);
                    yCompare = floor(c + (extendLength-1) * yStepSize);
                    
                    % Break if outside image
                    if xMeasure <= 0 || xMeasure > m ||...
                            yMeasure <= 0 || yMeasure > n
                        break
                    end
                    
                    % Continue if outside AF
                    if ~bw(xMeasure,yMeasure)
                        
                        extendSteps = extendSteps + 1;
                        
                        pixelDifference1 = ...
                            double(im1Blurred(xMeasure,yMeasure))-...
                            double(im1Blurred(xCompare,yCompare));
                        pixelDifference2 = ...
                            double(im2Blurred(xMeasure,yMeasure))-...
                            double(im2Blurred(xCompare,yCompare));
                        
                        if (pixelDifference1 <= 0) ||...
                                (extendSteps < minSteps) && ~done1
                            glow1(xMeasure,yMeasure) = 1;
                        else
                            done1 = 1;
                        end
                        
                        if (pixelDifference2 <= 0) ||...
                                (extendSteps < minSteps) && ~done2
                            glow2(xMeasure,yMeasure) = 1;
                        else
                            done2 = 1;
                        end
                        
                    end
                    
                    if (done1 && done2) || (extendSteps > maxSteps)
                        break
                    end
                    extendLength = extendLength+1;  
                end
                thMultiple = thMultiple + 1;
            end
        end
    end
end

im1GlowRemoved = im1;
im1GlowRemoved(bw==1) = 0;
im1GlowRemoved(glow1) = 0;

im2GlowRemoved = im2;
im2GlowRemoved(bw==1) = 0;
im2GlowRemoved(glow2) = 0;

