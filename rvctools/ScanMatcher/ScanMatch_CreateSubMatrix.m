function SubMatrix = ScanMatch_CreateSubMatrix(binX, binY, Threshold)
% SCANMATCH_CREATESUBMATRIX creates a substitution matrix based on the
% Euclidean distance between each RoI (bin). The Threshold value can be set
% by relating it to the variability of the saccade amplitude of your data 
% set and can be calculated as follow:
%
%	Threshold = 2 * std(saccade amplitude in pixels) / (Xres / Xbin)
%
% This means that the alignment algorithm will aim to align only regions 
% which are within the variability of the saccade amplitudes.
%
%   Part of the ScanMatch toolbox
%   Written by Filipe Cristino 
%   $Version: 1.00 $  $Date: 10/09/2009


% ---- Compute the Eucledian distance between each bin of the mask ----
mat = zeros(binX * binY, binX * binY);
indI = 1;
indJ = 1;
for i=1:binY % Need to find a better way to do this!
    for j=1:binX
        for ii=1:binY
            for jj=1:binX
               mat(indI,indJ) = sqrt((j-jj)^2 + (i-ii)^2);
               indI = indI+1;
            end
        end
        indI = 1;
        indJ = indJ+1;
    end
end

% ---- Inverse and Normalise the matrix with the Threshold
max_sub =  max(mat(:));
SubMatrix = abs(mat -max_sub) - (max_sub - Threshold);
end