function line = LsqLine(points)
x=points(1,:);
y=points(2,:);
n = length(x);

x_e = sum(x)/n;
y_e = sum(y)/n;

x_sum = sum(x);
y_sum = sum(y);

x2_sum = sum(x.^2);
y2_sum = sum(y.^2);

xy_sum = sum(x.*y);

alpha_x = x_sum^2 - y_sum^2 - n*x2_sum + n*y2_sum;
alpha_y = 2*x_sum*y_sum - 2*n*xy_sum;

alpha = atan2(alpha_y, alpha_x)/2;
r = x_e*cos(alpha)+y_e*sin(alpha);

if r<0 
    r = -r;
    if alpha < 0
        alpha = alpha + pi;
    elseif alpha > pi
        alpha = alpha - pi;
    end
end


line = [alpha;r];

end