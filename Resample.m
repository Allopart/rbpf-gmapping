%% Resampling function
function [x_new, q_new, idx] = Resample(x, q, resampling_strategy)

    Ns = length(q);  % Ns = number of particles

    switch resampling_strategy
       case 'multinomial_resampling'
          with_replacement = true;
          idx = randsample(1:Ns, Ns, with_replacement, q);

       case 'systematic_resampling'
          % this is performing latin hypercube sampling on q
          edges = min([0 cumsum(q)'],1); % protect against accumulated round-off
          edges(end) = 1;                 % get the upper edge exact
          u1 = rand/Ns;
          % this works like the inverse of the empirical distribution and returns
          % the interval where the sample is to be found
          [~, idx] = histc(u1:1/Ns:1, edges);
       otherwise
          error('Resampling strategy not implemented')
    end;

    x_new = x(:,idx);                    % extract new particles
    q_new = repmat(1/Ns, 1, Ns);          % now all particles have the same weight

end