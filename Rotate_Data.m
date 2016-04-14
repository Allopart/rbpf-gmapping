function data_out = Rotate_Data(data, pose)
    %% Rotates a set of data around one point and its angle
    
    theta=pose(3);
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    center = repmat([pose(1); pose(2)], 1, size(data, 2));
    so = R*data;                    % apply the rotation about the origin
    data_out = so + center;         % shift again so the origin goes back to the desired center of rotation
    
    
end