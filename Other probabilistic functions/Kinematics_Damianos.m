function [x_new] = Kinematics_Damianos(x,u,l)
 x_new = x + [cos(x(3)); sin(x(3)); tan(u(2))/l]*u(1);
end

