function [x_new] = Kinematics(x,u)
    %% Calculates the end position given an initial position and a movemnt command
    % Probabilistic Robotics -Thrun, Burgard, Fox 3rd Edition pg 127
    
    if u(2)==0
        u(2)=1e-10;     % If rotational value of command is exactly zero it gives errors
    end
    
    t_diff=1;                % Time step
    x_new = x + [-(u(1)/u(2))*sin(x(3))+(u(1)/u(2))*sin(x(3)+u(2)*t_diff); (u(1)/u(2))*cos(x(3))-(u(1)/u(2))*cos(x(3)+u(2)*t_diff); u(2)*t_diff];

    if x_new(3)<-2*pi
        x_new(3)=x_new(3)+2*pi;
    elseif x_new(3)>2*pi
        x_new(3)=x_new(3)-2*pi;
    end
end

