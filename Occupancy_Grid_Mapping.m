function [P_new, L_new] = Occupancy_Grid_Mapping(L_past, pose, scan, UsableArea, L0, gridsize, LRS) 
    %% Calculation of the occupancy grid map
    % Probabilistic Robotics-Thrun, Burgard, Fox 3rd Edition pg 286

    % Find area of interest
    x=round(pose(1));
    y=round(pose(2));
    L_new=L_past;
    
    for a=(x-UsableArea):(x+UsableArea)
         if a>=1 && a<=gridsize(1)
            for b=(y-UsableArea):(y+UsableArea)
                if (b)>=1 && (b)<=gridsize(2)
                    L_new(gridsize(2)-b, a)=L_past(gridsize(2)-b, a)+Inverse_Range_Sensor_Model(a,b,pose,scan, L0, UsableArea, LRS)- L0;
                end
            end
         end
    end
   
    % Calculate normalized probablity form log-odds
    P_new=1-(1./(1+exp(L_new)));
    

	%% Visulaize results: check if grid map is being built correctly (
%     figure(85)
%     set(gcf,'units','normalized','outerposition',[-1 0 1 1]);
%     clf;
%     hold on;
%     axis([0 gridsize(1) 0 gridsize(2)]);
%     imagesc(flipud(L_new))
%     scan_c=Robot2World(pose, Polar2Cart(scan));
%     plot(scan_c(1,:), scan_c(2,:), '.r');
%     hold off;
end

