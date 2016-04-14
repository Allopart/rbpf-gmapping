function [cart] = Polar2Cart(polar)
 %% Cnverts polar coordinats [alpha;r] to cartesian [x;y]

    alpha = polar(1,:);
    r = polar(2,:);
    x = zeros(1,size(polar,2));
    y = zeros(1,size(polar,2));

    for i = 1:size(polar,2)
    x(i) = r(i)*cos(alpha(i));
    y(i) = r(i)*sin(alpha(i));
    end

    cart = [x;y];
end

