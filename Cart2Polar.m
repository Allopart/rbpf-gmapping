function [polar] = Cart2Polar(cart)
    %% Cart2Polar converts cartesian coordinats [x;y] to polar [alpha;r]
    x = cart(1,:);
    y = cart(2,:);
    alpha = zeros(1,size(cart,2));
    r = zeros(1,size(cart,2));

    for i = 1:size(cart,2)
    alpha(i) = atan2(y(i),x(i));
    r(i) = sqrt(x(i)^2 + y(i)^2);
    end

    polar = [alpha;r];
end

