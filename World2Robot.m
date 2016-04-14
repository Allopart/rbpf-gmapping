function [x_robot] = World2Robot(x_world,x)
%ROBOT2WORLD Summary of this function goes here
%   Detailed explanation goes here
x_robot = zeros(size(x));
for i = 1:size(x,2)
x_robot(1,i) = cos(x_world(3))*(x(1,i)) + sin(x_world(3))*(x(2,i)) - x_world(1);
x_robot(2,i) = cos(x_world(3))*(x(2,i)) - sin(x_world(3))*(x(1,i)) - x_world(2);
end
end

