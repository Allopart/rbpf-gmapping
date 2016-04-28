function [xh,pf, P, L] = RBPF(x_true, y, LRS, u, a, pf, Nsamples, resampling_strategy)
    
    % x_true: true position of robot
    % y: True laser scan
    % LRS: laser scan inforamtion
    % u: Movement control
    % a: Odometry parameters found via Motion_Model_Velocity_Test
    % pf: pf class with all particle filter information
    % Nsamples: number of samples for the Particle Filter
    % resampling_startegy: DO NOT CHANGE    

    %% Variable Initialization

    k = pf.k;
    NPC = pf.NPC;
    q_past = pf.q(:,k-1);
    x_past= pf.particles(:,:,k-1);
    
    % Initialise  new variables
    x_initial_estimate=zeros(3, NPC);           % Initial estimate from previous Xt-1 and Ut
    yh=zeros(2, 361, NPC);                      % Laser scan estimate from x_initial_estimate
    x_new_sm=zeros(3, NPC);                     % New position after scan matching (ICP)
    x_sampled=zeros(3,Nsamples, NPC);           % Smaple position around every particle position
    p_xt=zeros(Nsamples, NPC);                  % Density function for odometry (Motion_Model_velocity)
    p_zt=zeros(Nsamples, NPC);                  % Density function for scans (Scan Matching)
    tau_xt=zeros(Nsamples, NPC);                % Optimal Proposal Distribution
    n=zeros(1,NPC);                             % Sample normalizer
    mu=zeros(3, NPC);                           % Mean for Gaussian distribution from all particle samples 
    cov=zeros(1, NPC);                          % Covariance for Gaussian distribution from all particle samples 
    q=zeros(size(q_past));                      % Particle weights
    x_final=zeros(3, NPC);                      % final estimate particle position after RBPF

    %% RBPF for every particle
    for i = 1:NPC

        %Estimate particle position from x_(t-1) after Ut
        x_initial_estimate(:,i) = Kinematics(x_past(:,i),u+sqrt(pf.R1)*randn(length(u),1));
        
        if  x_initial_estimate(3,i)>2*pi
            x_initial_estimate(3,i)=x_initial_estimate(3,i)-2*pi;
        elseif x_initial_estimate(3,i)<-2*pi
            x_initial_estimate(3,i)=x_initial_estimate(3,i)+2*pi;
        end
        
        %Scan Matcher returns estimate of laser given an initial particle position estimate
        yh(:,:,i) = Measurement_Estimate_From_Grid(x_initial_estimate(:,i),pf.P(:,:, i),LRS.Resolution,LRS.FoV, LRS.MaxDistance, pf.UsableArea);
             
        % Scan pre-processing prior to ICP algortithm
            y_help=y;
            yh_help=yh(:,:,i);
            
            % Use data inside Usable Area only
            y_help(y_help>pf.UsableArea)=nan;
            yh_help(yh_help>pf.UsableArea)=nan;
            
            % Chose data that whose angle is usable in both scans
            find_y_help=isnan(y_help(2,:));
            find_yh_help=isnan(yh_help(2,:));
            del_cols=or(find_y_help, find_yh_help);
            y_help(:,(del_cols(1,:)==1))=[];
            yh_help(:,(del_cols(1,:)==1))=[];
            
            % Transform to Polar coordinates and World frame
            y_help=Polar2Cart(y_help);
            yh_help=Polar2Cart(yh_help);
            y_help=Rotate_Data(y_help, x_true);
            yh_help=Robot2World([x_initial_estimate(1,i); x_initial_estimate(2,i); x_initial_estimate(3,i)], yh_help);
            
            % Add the particles position so that its also transformed  
            y_help=horzcat(y_help, x_true(1:2));
            yh_help=horzcat(yh_help, x_initial_estimate(1:2,i));
        

        failure=false;                               % Initialse ICP failure to FALSE
        
        try
        [TR,TT, data] = ICP(y_help, yh_help);        % ICP algortithm returns the translational and rotational matrices and the resulting data
        catch
            failure=true;
        end
        % Check if ICP failed
        if abs(TT(1))>50  || abs(TT(2))>50
            failure=true;
        end
        
        if failure
            disp('ICP failed')
            x_final(:,i)=x_initial_estimate(:,i);    % If ICP failed used last known position estimate
            q(i)=q_past(i);                          % do not change particle weights
        else
                     
            theta_new_sm=x_initial_estimate(3,i)+asin(TR(2,1));     % Modify theta angle according to the rotation given by ICP
            x_new_sm(:,i)=[data(1:2,end); theta_new_sm];            % Use rotated+translated initila_estimate as new (and better) positon
            
            %Plot grid-map, all positions and all scans (Black: True postion, Blue: Position and scan initial estimates, red:
            %Converged position and laser data
                figure(42)
                set(gcf,'units','normalized','outerposition',[1 0 1 1]);
                clf;
                caxis([0, 1]);
                colormap gray;
                colormap(flipud(colormap));
                hold on;
                imagesc(flipud(pf.P(:,:, i)));

                % Plot scans
                plot(y_help(1,:), y_help(2,:), '.g')
                plot(yh_help(1,:), yh_help(2,:), '.b')
                plot(data(1,:), data(2,:), '.r')

                % Plot positions
                plot(x_true(1), x_true(2), 'dk')
                plot(x_initial_estimate(1,i), x_initial_estimate(2, i), 'db')
                plot(x_new_sm(1,i), x_new_sm(2,i), 'dr')
                hold off;

            % Randomize samples around good estimate (x_new_sm) of particle position 
            theta_samples=x_new_sm(3,i);
            if theta_samples > 2*pi
                theta_samples=theta_samples-2*pi;
            elseif theta_samples < -2*pi
                theta_samples=theta_samples-2*pi;
            end    
            
            x_sampled(:,:, i) = [ normrnd([x_new_sm(1,i).*ones(1, Nsamples); x_new_sm(2,i).*ones(1, Nsamples)],0.5, [2 Nsamples]); theta_samples.*ones(1, Nsamples)];          
            
            for h=1:Nsamples
                              
                % Calculate odometry probability for that sample 
                p_xt(h,i)= Motion_Model_Velocity(x_sampled(:,h,i), u, x_past(:,i), a);
                              
                % Calculate scan matching probability for that sample 
                ytru=y;
                yest = Measurement_Estimate_From_Grid(x_sampled(:,h,i),pf.P(:,:, i),LRS.Resolution,LRS.FoV, LRS.MaxDistance, pf.UsableArea);
                yest(yest>pf.UsableArea)=nan;
                find_yest=isnan(yest(2,:));
                yest(:,(find_yest(1,:)==1))=[];
                ytru(:,(find_yest(1,:)==1))=[];
                p_zt(h,i) = Measurement_Model(ytru,yest,[1 1]);

                % Calculate optimal proposal distribution
                tau_xt(h,i)=abs(p_xt(h,i)*p_zt(h,i));

            end
            
            n(i)=sum(tau_xt(:,i));                                  % Calculate sample normalizer
            mu(1, i)=(1/n(i))*(x_sampled(1,:,i)*tau_xt(:,i));       % Calculate mu for gaussian distribution from samples
            mu(2, i)=(1/n(i))*(x_sampled(2,:,i)*tau_xt(:,i));
            mu(3, i)=theta_samples;                                 % Orientation of samples is constant (the one from ICP)

            % Calculate covariance for gaussian distribution from samples
            cov_sum=0;
            for h=1:Nsamples
                cov_sum=cov_sum+(x_sampled(1:2,h,i)-mu(1:2,i))'*(x_sampled(1:2,h,i)-mu(1:2,i))*tau_xt(h,i);
            end
            cov(i)=(1/n(i))*cov_sum;

            %Calculate new weights
            q(i)=q_past(i)*n(i);

            %Sample gaussian for final position of particle
            x_final(:,i)=[normrnd([mu(1,i); mu(2,i)],cov(i), [2 1]); mu(3, i)];
        end     
        
        %Update GridMap for taht particle
        [pf.P(:,:,i), pf.L(:,:,i)] = Occupancy_Grid_Mapping(pf.L(:,:,i), x_final(:,i), y, pf.UsableArea, pf.L0, pf.gridsize, LRS);
        
    end
    
    q=q/norm(q,1);                      % Normalize weights
    Neff = 1/sum(q.^2);                 % Calcuklate Neff for possible resampling
    
    
    %% Resampling step
    
    resample_percentage = 0.3;          % 0.5 is the best resample percentage proposed by Cyrill Stachniss
    Nt = resample_percentage*NPC;
    if Neff < Nt
        disp('Resampling');
        [x_final(:,:), q, idx] = Resample(x_final(:,:), q, resampling_strategy);
        P_save_resamp=pf.P(:,:,:);
        L_save_resamp=pf.L(:,:,:);
        for n=1:NPC
           pf.P(:,:,n)=P_save_resamp(:,:,idx(n));
           pf.L(:,:,n)=L_save_resamp(:,:,idx(n));
        end
    end
    
    pf.q(:,k) = q;                      % Save particle weights (might have changed during resampling)  
    
    %% Calculate estimate of true position (xh) given all weighted particles
    xh = zeros(3,1);
    for i=1:NPC
            xh = xh +(q(i).*x_final(:,i));
    end
    pf.xh(:, k)= xh;
    
    %% Store final position of all particles
    pf.particles(:,:,k) = x_final(:,:);
    
    % Make sure no errors have occured. If so probably becasue of p_xt or p_zt being equal to 0
    if isnan(x_final(1,i)) || isnan(x_final(2,i))
                1; %debugger
                x_final=x_initial_estimate;
    end
    
    
    %% Visualize results: Plot Particles 
%         figure(89)
%         set(gcf,'units','normalized','outerposition',[-1 0 1 1]);
%         clf;
%         hold on;
%         for g = 1:size(map,2)
%             plot([map(1,g) map(3,g)],[map(2,g) map(4,g)],'k')
%         end
%         plot(squeeze(pf.particles(1,:,pf.k-1))',squeeze(pf.particles(2,:,pf.k-1))','bo');
%         plot(squeeze(pf.particles(1,:,pf.k))',squeeze(pf.particles(2,:,pf.k))','k+');
%         plot(squeeze(pf.xh(1,pf.k))',squeeze(pf.xh(2,pf.k))','ro');
%         plot(x_true(1), x_true(2), 'go')
%         hold off;

end


