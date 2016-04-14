function p_xt= Motion_Model_Velocity_Test(u, x_p, a, gridsize)
%% Test to figure out which 'a' parameters suit you best
    % u: control command
    % x_p: initial position with orientation
    % a: [a_1 ... a_6 ] of tested parameters
    
    %the test will calculate the probability of the robot moving to every
    %possible cell in the grid generating a banana-shaped density function
    p_xt=zeros(gridsize(1), gridsize(2));
    for i=1:gridsize(1)
        for j=1:gridsize(2)-1
            p_xt(i,j) = Motion_Model_Velocity([j i x_p(3)], u, x_p, a) ; 
        end
    end
    
    % Visualize results
    figure(66)
    clf;
    hold on;
    plot(x_p(1), x_p(2), 'dr');
    surf(p_xt);
    view([0 -90]);
    hold off;

end