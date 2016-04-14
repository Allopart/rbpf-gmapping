function [q] = MeasurementModel(y,yh,sigma_v)
%MEASUREMENTMODEL Summary of this function goes here
%   Detailed explanation goes here

p_a = exp(-0.5*((y(1) - yh(1))/(sigma_v(1)))^2)/(sqrt(2*pi)*sigma_v(1));
p_r = exp(-0.5*((y(2) - yh(2))/(sigma_v(2)))^2)/(sqrt(2*pi)*sigma_v(2));
q = p_a*p_r;

end

