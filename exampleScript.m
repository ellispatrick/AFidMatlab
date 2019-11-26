%% Step 1: Read in images
ch1 = imread('1-1.tif');
ch2 = imread('1-2.tif');

%% Step 2: Obtain mask for analysis

% Read in mask
mask = logical(imread('mask.tif'));

% Generate mask
% sig = 2;
% ch1Blurred = imgaussfilt(ch1,sig);
% ch2Blurred = imgaussfilt(ch2,sig);
% 
% ch1Threshold = imbinarize(ch1Blurred,graythresh(ch1Blurred));
% ch2Threshold = imbinarize(ch2Blurred,graythresh(ch2Blurred));
%  
% mask = ch1Threshold & ch2Threshold;

% Size filter
minArea = 20;
maxArea = 10000;

mask = bwareafilt(mask,[minArea maxArea]);

%% Step 3: Autofluorescence classification

% Correlation only
% [maskAF, ch1AFRemoved, ch2AFRemoved] = afIdentifier(ch1, ch2, mask, 'Corr', 0.60);
% 
% Clustering with given k
% [maskAF, ch1AFRemoved, ch2AFRemoved] = afIdentifier(ch1, ch2, mask, 'k', 6);
% 
% Clustering with estimated k
[maskAF, ch1AFRemoved, ch2AFRemoved, kBest] = afIdentifier(ch1, ch2, mask,  'kAuto', 1, 'k', 20);
% 
% Clustering and correlation
% [maskAF, ch1AFRemoved, ch2AFRemoved] = afIdentifier(ch1, ch2, mask, 'k', 6, 'Corr', 0.60);
% OR
% [maskAF, ch1AFRemoved, ch2AFRemoved, kBest] = afIdentifier(ch1, ch2, mask,  'kAuto', 1, 'k', 20, 'Corr', 0.60);

%% Step 4: Glow Removal
[glow1, glow2, ch1GlowRemoved, ch2GlowRemoved] = glowIdentifier(ch1, ch2, maskAF, 'TraceSensitivity', 20, 'Sigma', 2);
