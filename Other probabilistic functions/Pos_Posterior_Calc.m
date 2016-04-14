
function p_xt = Pos_Posterior_Calc(x_p, u, gridsize, max_range)

p_xt=zeros(gridsize);

for i=x_p(1)-max_range:x_p(1)+max_range
    if i<=0 || i>gridsize(1)
        %skip
    else
        for j=x_p(2)-max_range:x_p(2)+max_range  
            if j<=0 || j>gridsize(2)
                %skip
            else
                x=[j i x_p(3)];
                p_xt(i,j) = Motion_Model_Velocity(x, u, x_p);
            end
    %        [p_xt1(i,j),p_xt2(i,j),p_xt3(i,j), mu(i,j)] = Motion_Model_Velocity(x, u_t, x_p,gridsize);
        end
    end
end

% figure(2)
% clf
% hold on
% surf(p_xt1);
% plot(x_p(1),x_p(2),'dr');
% hold off
% figure(3)
% clf
% hold on
% surf(p_xt2);
% plot(x_p(1),x_p(2),'dr');
% hold off
% 
% figure(4)
% clf
% hold on
% surf(p_xt3);
% plot(x_p(1),x_p(2),'dr');
% hold off
% 
% figure(5)
% clf
% hold on
% surf(p_xt3.*p_xt2.*p_xt1);
% plot(x_p(1),x_p(2),'dr');
% hold off

p_xt=p_xt./max(max(p_xt));

figure(6)
surf(p_xt);
hold on;
plot(x_p(1), x_p(2),'dr')
view(0,90)
hold off;
end