function p_xt = Motion_Model_Velocity(x, u, x_p, a) 
    %% Calculates the probability desity function of a specific position given an initial position and a movement command
    % Probabilistic Robotics-Thrun, Burgard, Fox 3rd Edition pg 123

    t_diff=1;                     % value of time steps
    x_1=x_p(1);
    y_1=x_p(2);
    theta_1=x_p(3);

    if theta_1>=pi
        theta_1=theta_1-2*pi;
    end
    if theta_1<=-pi
        theta_1=theta_1+2*pi;
    end

    x_2=x(1);
    y_2=x(2);
    theta_2=x(3);


    if x_2==22 && y_2==24
        1;
    end

    pos_cell_relative=[x_2-x_1, y_2-y_1];
    theta_cell=atan2(pos_cell_relative(2),pos_cell_relative(1));
    if theta_cell<theta_1 && (theta_cell+2*pi)>(theta_1+pi)
        h = Mirror_Point(pos_cell_relative,[tan(theta_1-pi/2) 0]);
        x_2=x_1+h(1);
        y_2=y_1+h(2);
    end

    v=u(1);
    w=rad2deg(u(2));

    mu= 0.5*(((x_1-x_2)*cos(theta_1)+(y_1-y_2)*sin(theta_1))/((y_1-y_2)*cos(theta_1)-(x_1-x_2)*sin(theta_1)));

    if mu==Inf
        mu=-abs(x_1-x_2);
    end
    if mu==-Inf
        mu=abs(x_1-x_2);
    end

    x_star=((x_1+x_2)/2)+mu*(y_1-y_2);
    y_star=((y_1+y_2)/2)+mu*(x_2-x_1);
    r_star=sqrt((x_1-x_star)^2+(y_1-y_star)^2);

    theta_diff1=atan2((y_2-y_star),(x_2-x_star));
    theta_diff2=atan2((y_1-y_star),(x_1-x_star));


    theta_diff=theta_diff1-theta_diff2;
    if theta_diff<-pi/2
        theta_diff=theta_diff+2*pi;
    end


    v_new=(theta_diff/t_diff)*r_star;
    w_new=theta_diff/t_diff;
    % theta_2=theta_1+w_new;
    gamma=((theta_2-theta_1)/t_diff)-w_new;

    p_xt1=Prob_Norm_Distribution(v-v_new,a(1)*v^2+a(2)*w^2);
    p_xt2=Prob_Norm_Distribution(w-w_new,a(3)*v^2+a(4)*w^2);
    p_xt3=Prob_Norm_Distribution(gamma,a(5)*v^2+a(6)*w^2);      % For simplification purposes we do not add a rotation at the end of the movement

    p_xt=p_xt1*p_xt2;%*p_xt3;

end