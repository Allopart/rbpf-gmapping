function scan_estimate = Measurement_Estimate_From_Grid(pose,grid_map,resol,FOV, max_range, UsableArea)
    %% Extracts  laser scan data from an occupancy grid map

    grid_map=flipud(grid_map);          % Map must be flipped because of the coordinates in which it is created
    object_limit = 0.817;               % Probability that says that an object is in that cell very certainly
	
	%Extract Pose
	x=pose(1);
	y=pose(2);
    
    if pose(3)<-2*pi
        pose(3)=pose(3)+2*pi;
    elseif pose(3)>2*pi
        pose(3)=pose(3)-2*pi;
    end
    
	theta = -pose(3); %% This negative is SUPER important
    
    
	%Extract Pose in grid coordinates
	xg = round(pose(1));
	yg = round(pose(2));
	gridsizex=size(grid_map,1);
	gridsizey=size(grid_map,2);

	% Number of scanning lines (deduced from scan width and resolution)
	num_of_scanning_lines=floor(FOV/resol)+1;
	resolrad=resol*pi/180.0;

    % initialise scan estimate to maximum scan range
	scan_estimate=[linspace(-resolrad*num_of_scanning_lines/2, resolrad*num_of_scanning_lines/2, num_of_scanning_lines); max_range*ones(1, num_of_scanning_lines)]; %Must make sure this matches with other scan arrays in diferent functions

    for a=(xg-UsableArea):(xg+UsableArea)
        if a>=1 && a<=gridsizex
            for b=(yg-UsableArea):(yg+UsableArea)
                if b>=1 && b<=gridsizey
                    if grid_map(b,a)>=object_limit

                        %Center of mass of cell
                        xc = a-0.5;
                        yc = b-0.5;

                        theta_cell=atan2(yc-y,xc-x);    % Angle form position cell to studied cell
                        if theta_cell<-pi
                            theta_cell=theta_cell+2*pi;
                        end
                        
                        theta_scan=theta_cell + theta;
                        if theta_scan<=-pi/2
                            theta_scan=theta_scan+2*pi;
                        end 

                        %Find estimate of number of scan closer to theta
                        rng=floor((theta_scan+resolrad*num_of_scanning_lines/2)/resolrad); %(max_angle+min_angle)/range=(1.5773+1.5773)/361

                        % Calculate distance
                        r=sqrt((xc-xg)^2+(yc-yg)^2);

                        %Calculates how many scans will go through a cell depending on how far it is
                        theta_rng=floor((atan(0.5/r)/resolrad)); 
                        
                        %Update that scan angle if its new distance is found to be samaller
                        for d=rng-theta_rng:rng+theta_rng
                            if d<=0 || d>=size(scan_estimate,2)
                                %skip
                            elseif r<scan_estimate(2,d)
                                scan_estimate(2,d)=r;
                            end
                        end		
                    end					
                end
            end
        end
    end   
    
    
    
    %% Visulaize results
%     figure(51)
%     clf;
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%     hold on;
%     imagesc(grid_map);
%     scan_e=Polar2Cart(scan_estimate(:,:));
%     scan_e=Robot2World([pose(1:2); pose(3)], scan_e);
%     plot(scan_e(1,:), scan_e(2,:), '.r');
%     plot(xg, yg, 'dr'); %mypos
%     quiver(xg, yg, 10*cos(theta), 10*sin(theta)),
%     hold off;


end