%SE2 Create planar translation and rotation transformation
%
% T = SE2(X, Y, THETA) is an SE(2) homogeneous transformation (3x3) 
% representing translation X and Y, and rotation THETA in the plane.
%
% T = SE2(XY) as above where XY=[X,Y] and rotation is zero
%
% T = SE2(XY, THETA) as above where XY=[X,Y]
%
% T = SE2(XYT) as above where XYT=[X,Y,THETA]
%
% See also TRANSL2, ROT2, ISHOMOG2, ISROT2, TRPLOT2.


% Copyright (C) 1993-2015, by Peter I. Corke
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

function t = se2(a, b, c, varargin)

    opt.deg = false;
    
    opt = tb_optparse(opt, varargin);
    
    if length(a) == 3
        x = a(1);
        y = a(2);
        th = a(3);
    elseif length(a) == 2
        x = a(1);
        y = a(2);
        if nargin < 2
            th = 0;
        else
            th = b;
        end
    else
        x = a;
        y = b;
        if nargin < 3
            th = 0;
        else
            th = c;
        end
    end

    if opt.deg 
        th = th * pi/180.0;
    end
    cth = cos(th);
    sth = sin(th);
    R = [cth -sth; sth cth];

    t = [R [x; y]; 0 0 1];
