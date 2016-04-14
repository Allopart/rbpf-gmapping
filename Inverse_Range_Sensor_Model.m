function result = Inverse_Range_Sensor_Model(ii,jj,pose,scan, L0, max_range, LRS)
    %% Calculation of the inverse model returning the likelihood of occupancy
    % Probabilistic Robotics -Thrun, Burgard, Fox 3rd Edition pg 288 


    alpha= 2;               % Thickness of obsatcles
    beta= 0.1;              % Width of sensor beam

    %Center of mass of cell
    xc = ii-0.5;
    yc = jj-0.5;        
    
    if pose(3)>2*pi
        pose(3)=pose(3)-2*pi;
    elseif pose(3)<-2*pi
        pose(3)=pose(3)+2*pi;
    end
    
    r=sqrt((xc-pose(1))^2+(yc-pose(2))^2);          % Calculate distance between psotion cell and studied cell
    theta=atan2(yc-pose(2), xc-pose(1))-pose(3);    % Calculate angle between psotion cell and studied cell
    
    if theta>2*pi
        theta=theta-2*pi;
    elseif theta<-2*pi
        theta=theta+2*pi;
    end
    
    if theta<-pi/2
        theta=theta+2*pi;
    end

    %Find estimate of number of scan closer to theta
    rng=floor((theta+LRS.MaxAngle)/(2*LRS.MaxAngle/(LRS.FoV/LRS.Resolution +1))); %(max_angle+min_angle)/range=(1.5773+1.5773)/361
    
    if rng>=size(scan,2)-1
        rng=size(scan,2)-1;
    end
    if rng<2
        rng=2;
    end

    %Search around the number of scan
    %Argmin part
    argmin=abs(theta-scan(1,1));
    kmin=1;

    for d=rng-1:rng+1
       a = abs(theta-scan(1,d));
       if a<argmin
           argmin=a;
           kmin=d;
       end
    end


    if (r > min(max_range, scan(2,kmin)+alpha/2)) || (abs(theta-scan(1,kmin)) > beta/2)
        result=L0;
    elseif (scan(2,kmin)< max_range) && (abs(r-scan(2,kmin))<alpha/2)
        result=1;
    elseif (r<=scan(2,kmin))
        result=-1;
    else
        error('Grid occupancy failed');
    end
  
    
    
    

end