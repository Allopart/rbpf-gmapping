function q = Likelihood_Field_Range_Finder_Model(scan, pose, map, sigma_v, max_range)
	
    z_max=max_range;
    q=1;
    
    x=pose(1);
    y=pose(2);
    theta=pose(3);
    
    x_sens=0;
    y_sens=0;
    theta_sens=0;
    
    x_sm=zeros(1, size(scan, 2));
    y_sm=zeros(1, size(scan, 2));
    
    p_random=[(1/z_max)*ones(1,z_max-1) 0];
    


%     scan_cart=Polar2Cart(scan);

    for k=1:size(scan, 2)
        if scan(2,k)~= z_max
            x_sm(k)=x+x_sens*cos(theta)-y_sens*sin(theta)+scan(2,k)*cos(theta+theta_sens);
            y_sm(k)=y+y_sens*sin(theta)+x_sens*sin(theta)+scan(2,k)*sin(theta+theta_sens);
            x_sm(k)=round(x_sm(k));
            y_sm(k)=round(y_sm(k));

            d_min=z_max;
            for s=(x_sm(k)-z_max):(x_sm(k)+z_max)
                if s>0 && s<size(map, 1)
                    for t=(y_sm(k)-z_max):(y_sm(k)+z_max)
                        if t>0 && t<size(map, 2)
                            if map(s,t)>0.75
                                d=abs((x_sm(k)-s)^2+(y_sm(k)-t)^2);  %Maybe x and y need to be switched
                                if d<d_min
                                    d_min=d;
                                end
                            end
                        end
                    end
                end
            end
            
            z_hit=0.6;
            z_random=0.2;
            g_max=0.2;
            
            q=q*(normpdf(d_min, sigma_v(1)));
        end
    end
end

%% 
% z=50;
% z_max=80;
% p_max=[zeros(1,z_max-1) 1];
% p_hit=pdf('Normal', (0:1:z_max-1), z, sigma_v(2));
% p_random=[(1/z_max)*ones(1,z_max-1) 0];
% 
% p_tot=1/3*p_max + 1/3*p_hit + 1/3*p_random;
% plot(p_tot);