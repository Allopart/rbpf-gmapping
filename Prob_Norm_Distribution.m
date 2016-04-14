function p = Prob_Norm_Distribution(a,b)
    %% Calculates the probability distribution
    
    p=(1/sqrt(2*pi*b))*exp(-0.5*(a^2/b));
end