%DOCKFIGS Control figure docking in the GUI
%
% dockfigs causes all new figures to be docked into the GUI
%
% dockfigs(1) as above.
%
% dockfigs(0) causes all new figures to be undocked from the GUI


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

function dockfigs(arg)
    if nargin == 0
        arg = 1;
    end
    
    if arg
        set(0, 'DefaultFigureWindowStyle', 'docked')
    else
        set(0, 'DefaultFigureWindowStyle', 'normal')
    end
