function q = Measurement_Model(ytru, yest, sigma_v)
%% Calculates the probability of one point in the real scan matching to that same point in the estimated scan
q=1;
for i=1:size(ytru(2,:))
    p_a = exp(-0.5*((ytru(1,i) - yest(1,i))/(sigma_v(1)))^2)/(sqrt(2*pi)*sigma_v(1));
    p_r = exp(-0.5*((ytru(2,i) - yest(2,i))/(sigma_v(2)))^2)/(sqrt(2*pi)*sigma_v(2));
    p = p_a*p_r;
    q=q*p;
end