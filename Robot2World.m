function [x_world] = Robot2World(x_robot,x)
    %% Change of coordinates from robot fram to world frame

    x_world = zeros(size(x));
    for i = 1:size(x,2)
    x_world(1,i) = cos(x_robot(3))*(x(1,i)) - sin(x_robot(3))*(x(2,i)) + x_robot(1);
    x_world(2,i) = cos(x_robot(3))*(x(2,i)) + sin(x_robot(3))*(x(1,i)) + x_robot(2);
    end
end

