function x_mirror= Mirror_Point(x, line)
    %% Mirrors a certain point over a line 
    
    d=(x(1)+(x(2)-line(2))*line(1))/(1+line(1)^2);
    x_mirror(1)=2*d-x(1);
    x_mirror(2)=2*d*line(1)-x(2)+2*line(2);
end