%ISROT2 Test if SO(2) rotation matrix
%
% ISROT2(R) is true (1) if the argument is of dimension 2x2 or 2x2xN, else false (0).
%
% ISROT2(R, 'valid') as above, but also checks the validity of the rotation
% matrix.
%
% Notes::
% - A valid rotation matrix has determinant of 1.
%
% See also ISHOMOG2, ISVEC.



% Copyright (C) 1993-2014, by Peter I. Corke
%
% This file is part of The Robotics Toolbox for MATLAB (RTB).
% 
% RTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% RTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with RTB.  If not, see <http://www.gnu.org/licenses/>.
%
% http://www.petercorke.com

function h = isrot2(r, dtest)

    d = size(r);
    if ndims(r) >= 2
        h =  all(d(1:2) == [2 2]);

        if h && nargin > 1
            h = abs(det(r) - 1) < eps;
        end
    else
        h = false;
    end
