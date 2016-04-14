function p_zt = Scan_Matching( true_scan,pose,map,UsableArea)
    %% Calculates the probability function for a specific location given it own map and the true scan that you should see
    
    map=flipud(map);
    
	%Extract Pose in grid coordinates
	xg = round(pose(1));
	yg = round(pose(2));
    theta = pose(3);
	gridsizex=size(map,1);
	gridsizey=size(map,2);
    
    theta_m=[theta*ones(1, size(true_scan, 2)); zeros(1, size(true_scan, 2))];
    grid_pos= [xg*ones(1, size(true_scan, 2)) ; yg*ones(1, size(true_scan, 2))];

    p_cell=zeros(1,size(true_scan, 2));
    
    %Find indexes of scans that are close to max_range
    find=true_scan(2,:)>UsableArea;
    true_scan(:,(find==1))=[];
    pos_scan=Polar2Cart(true_scan + repmat([-theta; 0], 1, size(true_scan, 2))) + repmat([xg; yg], 1, size(true_scan, 2));
    pos_scan(pos_scan<=1)=1;
    pos_scan(pos_scan(1,:)>gridsizey)=gridsizey;
    pos_scan(pos_scan(2,:)>gridsizex)=gridsizex;
    pos_scan=round(pos_scan);
    
    for i=1:size(pos_scan,2)
        p_cell(i)=map(pos_scan(1,i), pos_scan(2,i));   
    end
    
    p_zt=sum(p_cell, 'omitnan');%/count;
%     p_zt= prod(p_cell);
    
    
%% Visualize results: 
    figure(9)
    hold on;
    imagesc(map);
    plot(pos_scan(1,:), pos_scan(2,:), '.r');
    plot(xg, yg, 'dk');
    hold off;
    
end