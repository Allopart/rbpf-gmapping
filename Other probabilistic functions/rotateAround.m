function output=rotateAround(image, pointY, pointX, angle, varargin)
% ROTATEAROUND rotates an image.
%   ROTATED=ROTATEAROUND(IMAGE, POINTY, POINTX, ANGLE) rotates IMAGE around
%   the point [POINTY, POINTX] by ANGLE degrees. To rotate the image
%   clockwise, specify a negative value for ANGLE.
%
%   ROTATED=ROTATEAROUND(IMAGE, POINTY, POINTX, ANGLE, METHOD) rotates the
%   image with specified method:
%       'nearest'       Nearest-neighbor interpolation
%       'bilinear'      Bilinear interpolation
%       'bicubic'       Bicubic interpolation
%    The default is fast 'nearest'. Switch to 'bicubic' for nicer results.
%
%   Example
%   -------
%       imshow(rotateAround(imread('eight.tif'), 1, 1, 10));
%
%   See also IMROTATE, PADARRAY.

%   Contributed by Jan Motl (jan@motl.us)
%   $Revision: 1.2 $  $Date: 2014/05/01 12:08:01 $

% Parameter checking.
numvarargs = length(varargin);
if numvarargs > 1
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 1 optional input');
end
optargs = {'nearest'};    % Set defaults for optional inputs
optargs(1:numvarargs) = varargin;
[method] = optargs{:};    % Place optional args in memorable variable names

% Initialization.
[imageHeight imageWidth ~] = size(image);
centerX = floor(imageWidth/2+1);
centerY = floor(imageHeight/2+1);

dy = centerY-pointY;
dx = centerX-pointX;

% How much would the "rotate around" point shift if the 
% image was rotated about the image center. 
[theta, rho] = cart2pol(-dx,dy);
[newX, newY] = pol2cart(theta+angle*(pi/180), rho);
shiftX = round(pointX-(centerX+newX));
shiftY = round(pointY-(centerY-newY));

% Pad the image to preserve the whole image during the rotation.
padX = abs(shiftX);
padY = abs(shiftY);

padded = padarray(image, [padY padX]);

% Rotate the image around the center.
rot = imrotate(padded, angle, method, 'crop');

% Crop the image.
output = rot(padY+1-shiftY:end-padY-shiftY, padX+1-shiftX:end-padX-shiftX, :);

