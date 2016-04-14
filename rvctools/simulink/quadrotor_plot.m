

% Copyright (C) 1993-2014, by Peter I. Corke
%
% This file is part of The Robotics Toolbox for Matlab (RTB).
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
function [sys,x0,str,ts] = quadrotor_plot(t,x,u,flag,s,plot,enable,vehicle)
    % Flyer plot, lovingly coded by Paul Pounds, first coded 17/4/02
    % version 2 2004 added scaling and ground display
    % version 3 2010 improved rotor rendering and fixed mirroring bug
    %
    % Displays X-4 flyer position and attitude in a 3D plot.
    % GREEN ROTOR POINTS NORTH
    % BLUE ROTOR POINTS EAST
    
    % PARAMETERS
    % s defines the plot size in meters
    % swi controls flyer attitude plot; 1 = on, otherwise off.
    
    % INPUTS
    % 1 Center X position
    % 2 Center Y position
    % 3 Center Z position
    % 4 Yaw angle in rad
    % 5 Pitch angle in rad
    % 6 Roll angle in rad
    
    % OUTPUTS
    %   None
    ts = [-1 0];
    
    if ~isfield(vehicle, 'nrotors')
        vehicle.nrotors = 4;    % sensible default for quadrotor function
    end
    
    switch flag,
        case 0
            [sys,x0,str,ts] = mdlInitializeSizes(ts,plot,enable); % Initialization
        case 3
            sys = mdlOutputs(t,u,s,plot,enable, vehicle); % Calculate outputs
        case {1,2, 4, 9} % Unused flags
            sys = [];
        otherwise
            error(['unhandled flag = ',num2str(flag)]); % Error handling
    end
    
    
    % Initialize
function [sys,x0,str,ts] = mdlInitializeSizes(ts,plot,enable)
    % Call simsizes for a sizes structure, fill it in, and convert it
    % to a sizes array.
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 0;
    sizes.NumInputs      = 6;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [];
    str = [];          % Set str to an empty matrix.
    ts = [0.05 0];
    
    if enable == 1
        figure(plot);
        clf;
        %colordef(1,'none');
    end
    % End of mdlInitializeSizes.
    
    
function sys = mdlOutputs(t,u,s, plot, enable, quad)
    global a1s b1s
    
    % not quite sure what this is about -- PIC
    if numel(a1s) == [0];
        a1s = zeros(1, quad.nrotors);
        b1s = zeros(1, quad.nrotors);
    end
    
    % vehicle dimensons
    d = quad.d; %Hub displacement from COG
    r = quad.r; %Rotor radius

    for i = 1:quad.nrotors
        theta = (i-1)/quad.nrotors*2*pi;
        %   Di      Rotor hub displacements (1x3)
        % first rotor is on the x-axis, clockwise order looking down from above
        D(:,i) = [ d*cos(theta); d*sin(theta); 0];
        scal = s(1)/4;
        %Attitude center displacements
        C(:,i) = [ scal*cos(theta); scal*sin(theta); 0];
    end
    
    if enable == 1
        %draw ground
        figure(plot);
        clf;
        if length(s) == 1
            axis([-s s -s s 0 s]);
        else
            axis([-s(1) s(1) -s(1) s(1) 0 s(2)])
            s = s(1);
        end
        hold on;
        
        % plot the ground boundaries and the big cross
        plot3([-s -s],[s -s],[0 0],'-b')
        plot3([-s s],[s s],[0 0],'-b')
        plot3([s -s],[-s -s],[0 0],'-b')
        plot3([s s],[s -s],[0 0],'-b')
        plot3([s -s],[-s s],[0 0],'-b')
        plot3([-s s],[-s s],[0 0],'-b')
        
        %READ STATE
        z = [u(1);u(2);u(3)];
        n = [u(4);u(5);u(6)];
        
        %PREPROCESS ROTATION MATRIX
        phi = n(1);    %Euler angles
        the = n(2);
        psi = n(3);
        
        R = [cos(the)*cos(phi) sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi) cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi);   %BBF > Inertial rotation matrix
            cos(the)*sin(phi) sin(psi)*sin(the)*sin(phi)+cos(psi)*cos(phi) cos(psi)*sin(the)*sin(phi)-sin(psi)*cos(phi);
            -sin(the)         sin(psi)*cos(the)                            cos(psi)*cos(the)];
        
        %Manual Construction
        %Q3 = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];   %Rotation mappings
        %Q2 = [cos(the) 0 sin(the);0 1 0;-sin(the) 0 cos(the)];
        %Q1 = [1 0 0;0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
        %R = Q3*Q2*Q1;    %Rotation matrix
        
        %CALCULATE FLYER TIP POSITONS USING COORDINATE FRAME ROTATION
        F = [1 0 0;0 -1 0;0 0 -1];
        
        %Draw flyer rotors
        t = [0:pi/8:2*pi];
        for j = 1:length(t)
            circle(:,j) = [r*sin(t(j));r*cos(t(j));0];
        end
        
        for i = 1:quad.nrotors
            hub(:,i) = F*(z + R*D(:,i)); %points in the inertial frame
            
            q = 1; %Flapping angle scaling for output display - makes it easier to see what flapping is occurring
            Rr = [cos(q*a1s(i))  sin(q*b1s(i))*sin(q*a1s(i)) cos(q*b1s(i))*sin(q*a1s(i));   %Rotor > Plot frame
                0              cos(q*b1s(i))               -sin(q*b1s(i));
                -sin(q*a1s(i)) sin(q*b1s(i))*cos(q*a1s(i)) cos(q*b1s(i))*cos(q*a1s(i))];
            
            tippath(:,:,i) = F*R*Rr*circle;
            plot3([hub(1,i)+tippath(1,:,i)],[hub(2,i)+tippath(2,:,i)],[hub(3,i)+tippath(3,:,i)],'b-')
        end
        
        %Draw flyer
        hub0 = F*z;  % centre of vehicle
        for i = 1:quad.nrotors
            % line from hub to centre plot3([hub(1,N) hub(1,S)],[hub(2,N) hub(2,S)],[hub(3,N) hub(3,S)],'-b')
            plot3([hub(1,i) hub0(1)],[hub(2,i) hub0(2)],[hub(3,i) hub0(3)],'-b')
            
            % plot a circle at the hub itself
            plot3([hub(1,i)],[hub(2,i)],[hub(3,i)],'o')
        end
        
        % plot the vehicle's centroid on the ground plane
        plot3([z(1) 0],[-z(2) 0],[0 0],'--k')
        plot3([z(1)],[-z(2)],[0],'xk')

        % label the axes
        xlabel('x');
        ylabel('y');
        zlabel('z (height above ground)');
    end
        
    sys = [];
    % End of mdlOutputs.
