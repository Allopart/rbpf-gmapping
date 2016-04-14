function mask = ScanMatch_GridMask(X,Y, binX, binY)
% SCANMATCH_GridMask create a grid RoI mask with ROI
% values going from 1 to binX * binY.
%   (X,Y) are the size of the mask in pixels
%   (binX, binY) are the number of bins
%
%   Part of the ScanMatch toolbox
%   Written by Filipe Cristino 
%   $Version: 1.00 $  $Date: 10/09/2009

% Create the RoI mask 
array = reshape(1:binX*binY, binX, binY)';

% resize the mask to be at the right output size without using 'imresize'
m = double(binX) / X;
n = double(binY) / Y;
[Xi Yi] = meshgrid(1:binX, 1:binY);
[Xn Yn] = meshgrid(floor(1:m:binX+1-m), floor(1:n:binY+1-n));
mask = interp2(Xi, Yi, array,Xn, Yn,'nearest');

end