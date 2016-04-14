%SerialLink.JACOBN Jacobian in end-effector frame
%
% JN = R.jacobn(Q, options) is a 6xN Jacobian matrix for the robot in 
% pose Q. The manipulator Jacobian matrix maps joint velocity to 
% end-effector spatial velocity V = J0*QD in the end-effector frame.
%
% Options::
% 'trans'   Return translational submatrix of Jacobian
% 'rot'     Return rotational submatrix of Jacobian 
%
% Notes::
% - this Jacobian is often referred to as the geometric Jacobian
%
% Reference::
% 	Paul, Shimano, Mayer,
% 	Differential Kinematic Control Equations for Simple Manipulators,
% 	IEEE SMC 11(6) 1981,
% 	pp. 456-460
%
% See also SerialLink.jacob0, delta2tr, tr2delta.

% Ryan Steindl based on Robotics Toolbox for MATLAB (v6 and v9)
%
% Copyright (C) 1993-2011, by Peter I. Corke
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

function J = jacobn(robot, q, varargin)

    opt.trans = false;
    opt.rot = false;
    
    opt = tb_optparse(opt, varargin);
    
	n = robot.n;
	L = robot.links;		% get the links

	J = [];
	U = robot.tool;
	for j=n:-1:1,
		if robot.mdh == 0,
			% standard DH convention
			U = L(j).A(q(j)) * U;
		end
		if L(j).RP == 'R',
			% revolute axis
			d = [	-U(1,1)*U(2,4)+U(2,1)*U(1,4)
				-U(1,2)*U(2,4)+U(2,2)*U(1,4)
				-U(1,3)*U(2,4)+U(2,3)*U(1,4)];
			delta = U(3,1:3)';	% nz oz az
		else
			% prismatic axis
			d = U(3,1:3)';		% nz oz az
			delta = zeros(3,1);	%  0  0  0
		end
		J = [[d; delta] J];

		if robot.mdh ~= 0,
			% modified DH convention
			U = L(j).A(q(j)) * U;
		end
    end
    
    if opt.trans
        J0 = J0(1:3,:);
    elseif opt.rot
        J0 = J0(4:6,:);
    end
