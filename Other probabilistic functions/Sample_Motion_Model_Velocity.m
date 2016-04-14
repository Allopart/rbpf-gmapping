function [x_new] = Sample_Motion_Model_Velocity(x_p,u, a, l)
    v=u(1);
    w=deg2rad(u(2));
    t=1;
    
        
    v_new=v+Sample_Normal_Distribution(a(1)*v^2+a(2)*w^2);
    w_new=w+Sample_Normal_Distribution(a(3)*v^2+a(4)*w^2);
    gamma=Sample_Normal_Distribution(a(5)*v^2+a(6)*w^2);
    
    x_new(1)=x_p(1)-(v_new/w_new)*sin(x_p(3))+(v_new/w_new)*sin(x_p(3)+w_new*t);
    x_new(2)=x_p(2)+(v_new/w_new)*cos(x_p(3))-(v_new/w_new)*cos(x_p(3)+w_new*t);
    x_new(3)=x_p(3)+w_new*t+gamma*t;
    x_new(3) = mod(x_new(3) + pi,2*pi)-pi; %[-pi,pi]
end