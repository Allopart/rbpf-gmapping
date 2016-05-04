    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% RAO BLACKWELLIZED PARTICLE FILTER FOR GRID-BASED FAST SALM SIMULATOR
% 
% Based on 'Probabilistic Robotics' by Thrun, Burgard and Fox and several
% papers by Cyril Stachniss
% 
% Adrian Llopart Maurin 15.04.16
% PhD in the AUT Group at the Technical University of Denmark
% adllo@elektro.dtu.dk
% 
% 
% DESCRIPTION
% Simulation of a robot moving through a determined path in a predefined map and the 
% resulting occupancy grid of a particle filter which only takes the control commands
% and true laser scans as inputs
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clean up
clear variables
close all
clc

addpath(genpath('rvctools'))
rng(1)                                          % Adjust the noise realization here

%% Simulation Parameters

% General parameters
Ts = 1;                                         % Sampling time
n_cell=10;                                      % Number of cells per meter
NPC = 3;                                        % Number of particles used
usable_area = 4*n_cell;                         % Use laser data in 4 meter radius
sigma_v = [0.01 0.01];                          % Standard deviation of sensor (Holds for up to 20m) (m)
R1 = [0.1 0;0 0.001];                           % Variance matrix for odometry 
Nsamples = 10;                                  % Number of samples for Particle Filter
max_speed=4;                                    % max_speen*n_cell cm/timestep

% LRS_Sensor Structure
LRS.MaxDistance         = 80*n_cell;            % Maximum sensor measurement distance (m/10)
LRS.Resolution          = 0.5;                  % Sensor resolution (degrees)
LRS.FoV                 = 180;                  % Sensor field of view (degrees)
LRS.MaxAngle            = 1.5773;               % extreme angle of laser scanner


%% Map Generation

%Consists of (x,y) start and end points that define a line in the space
map = zeros(4,24);
%Wall Exterior
map(:,1) = [80; 20; 210; 20];
map(:,2) = [210; 20; 210; 230];
map(:,3) = [210; 230; 80; 230];
map(:,4) = [80; 230; 80; 170];
map(:,5) = [80; 170; 20; 170];
map(:,6) = [20; 170; 20; 80];
map(:,7) = [20; 80; 80; 80];
map(:,8) = [80; 80; 80; 20];
%Wall Interior
map(:,9) = [110; 50; 180; 50];
map(:,10) = [180; 50; 180; 200];
map(:,11) = [180; 200; 110; 200];
map(:,12) = [110; 200; 110; 170];
map(:,13) = [110; 170; 170; 170];
map(:,14) = [170; 170; 170; 80];
map(:,15) = [170; 80; 110; 80];
map(:,16) = [110; 80; 110; 50];
%Wall Interior
map(:,17) = [50; 110; 140; 110];
map(:,18) = [140; 110; 140; 140];
map(:,19) = [140; 140; 50; 140];
map(:,20) = [50; 140; 50; 110];


% % Plot map
% figure(1)
% clf
% hold on
% for d = 1:size(map,2)
%     plot([map(1,d) map(3,d)],[map(2,d) map(4,d)],'k')
% end
% 
% % First Path Generation
% via = [105,35; 185,35; 185,40; 190,45; 195,50; 195,203; 190,207; 185,210; 180,215; 105,215; 103,195; 100,190; 95,185; ...
%      95,170; 100,165; 105,160; 110,155; 140,155; 145,150; 150,145; 155,140; 155,110; 150,105; 145,100; 140,95; 50,95; ...
%      45,100: 40,105; 35,110; 35,140; 40,145; 45,150; 50,155; 110,155; 140,155; 145,150; 150,145; 155,140; 155,110; ...
%      150,105; 145,100; 140,95; 110,95; 105,90; 100,85; 95,80;95,35];
% %Controller Parameters
% d = 1; 
% Kp = 0.15;
% Kh = 0.5;

%% Second path generation 
via = [100,35; 185,35; 185,40; 190,45; 195,50; 195,200; 190,205; 185,210; 180,215; 110,215; 105,195; 100,190; 95,185; ...
    95,170; 90,165; 85,160; 80, 155; 50,155; 45, 150; 40,145; 35,140; 35,110; 40,105; 45, 100; 50,95; 140,95; 145,100; ...
    150,105; 155,110; 155,140; 150,145; 145,150; 140,155; 50,155; 45, 150; 30,145; 35,140; 35,110; 40,105; ...
    45,100; 50,95; 80,95; 85,90; 90,85; 95,80; 95,35]; 

%Controller Parameters
d = 1; 
Kp = 0.4;
Kh = 0.9;

%% Generation of a high order polynomial with the drive-through
% points and a desired velocity as boundary conditions.

path = mstraj(via,[max_speed, max_speed],[],[via(1,1) via(1,2)],Ts,0); 

switch max_speed                                %Calculated Model Velocity parameters (a) via Motion_Model_Velocity_Test
    case 2
        a= [0.05 0.01 0.01 0.05 0.0001 0.0001];
    case 3
        a= [0.1 0.001 0.001 0.1 0.0001 0.0001];
    case 4
        a=[0.1 0.01 0.01 0.1 0.0001 0.0001]; 
    case 5
        a=[0.1 0.001 0.001 0.1 0.0001 0.0001]; 
end       


%% Grid occupancy initialisation
Nstep=size(path,1);
security=1.2;                                   % Initialise the grid 120% bigger than the map actually is
xg=max(max(map))*security;                      % Length
yg=xg;                                          % Width

L=0.5*ones(yg,xg);                              % Initialise the probability of grid map

%% Data Initialization

%Allocating sufficient space for the variables.
Nsim = size(path,1);
u = zeros(2,Nsim);
x = zeros(3,Nsim);

x0 = [via(1,1),via(1,2),0];                     % Initialize robot pose at the start of the path.

%Particle Filter Structure
 pf.k               = 1;                        % Iteration number
 pf.NPC             = NPC;                      % Number of particles
 pf.q               = zeros(pf.NPC, Nsim);      % Weights
 pf.particles       = zeros(3, pf.NPC, Nsim);   % Particles position
 pf.sigma_v         = sigma_v;                  % Propagating the noise of the sensor
 pf.P               = zeros(yg, xg, pf.NPC);    % Actual probability of grid map
 pf.L               = zeros(yg, xg, pf.NPC);    % Log odds probability 
 pf.L0              = 0;                        % L0 value for inverse_range_sensor_model
 pf.gridsize        = [xg yg];                  % Gridsize
 pf.xh              = zeros(3, Nsim);           % Estiamte of true position
 pf.R1              = R1;                       % Variance matrix for odometry 
 pf.UsableArea      = usable_area;              % Laser scanner usable area

%% Time step = 1 initialization

tic
h = timebar('Progress','FastSLAM Simulation');  % Start timebar

x(:,1)=x0;                                      % Initialise particle position
pf.xh(:,1)=x0;                                  % Initialise estiamte
y = LaserScanNoise(pf.xh(:,1),map,LRS,sigma_v); % Obtain true laser data from map

% First laser scan data is asserted as occupancy grid TWICE for a very good initialisatio
[P_initial, L_initial] = Occupancy_Grid_Mapping(L, pf.xh(:,1), y, pf.UsableArea, pf.L0, [xg yg], LRS);
[P_initial, L_initial] = Occupancy_Grid_Mapping(L_initial, pf.xh(:,1), y, pf.UsableArea, pf.L0, [xg yg], LRS);

% Initialise all variables at timestep = 1
for i=1:NPC
pf.P(:,:,i)=P_initial;
pf.L(:,:,i)=L_initial;
pf.q(i,1)=1/NPC;                                % Weights are initialised uniformly
pf.particles(:, i, 1)=x0;                       % All particles start from the same position
end



%% Start simulation at second timestep

for k = 2:Nstep
    
    timebar(h, k/Nstep)

    % True System
        % Compute Error Signals
            e = sqrt((path(k,1)-x(1,k-1))^2 + (path(k,2)-x(2,k-1))^2) - d;
            th = atan2((path(k,2)-x(2,k-1)),(path(k,1)-x(1,k-1)));
        % Control Signal
            u(1,k-1) = Kp*e;
            u(2,k-1) = Kh*(angdiff(th,x(3,k-1)));
        % Kinematics
            x(:,k) = Kinematics(x(:,k-1),u(:,k-1));
        % True LRS data
            y = LaserScanNoise(x(:,k),map,LRS,sigma_v);         
            
     % The Rao-blackwellized Particle Filter for Grid Based fastSlam
            pf.k = k;
            [pf.xh(:,k),pf] = RBPF(x(:,k),y,LRS,u(:,k-1),a,pf,Nsamples,'multinomial_resampling');
        
     % Visulaization of true position moving through map
%         [pf.P(:,:,1), pf.L(:,:,1)] = Occupancy_Grid_Mapping(pf.L(:,:,1), x(:,k), y, pf.UsableArea, pf.L0, [xg yg]);
%         figure(3)
%         set(gcf,'units','normalized','outerposition',[1 0 1 1]);
%         clf    
%         imagesc(flipud(pf.P(:,:,1)));
%         caxis([0, 1]);
%         colormap gray;
%         colormap(flipud(colormap));
%         hold on;
%         grid on;
%         plot(x(1,k),x(2,k), 'dr');
%         hold off;
    
end
toc
close(h)


%% Grid end result visualization for all particles: They should be very similar due to resampling stages
figure(88)
for particles=1:pf.NPC   
    subplot(2,3,particles)       
    imagesc(pf.P(:,:,particles));
    caxis([0, 1]);
    colormap gray;
    colormap(flipud(colormap));
    grid on;
    str = sprintf('Grid map for particle %d',particles);
    title(str);
    hold off;
end

%% State Plots
figure(2)
subplot(3,1,1)
hold on
h1 = plot(squeeze(pf.particles(1,:,:))','y');
h2 = plot(x(1,:));
h3 = plot(pf.xh(1,:),'k');
grid on
title('X Estimation')
xlabel('Samples')
ylabel('Position (dm)')
legend([h1(1) h2 h3],'Particles','X','X_h')

subplot(3,1,2)
hold on
h1 = plot(squeeze(pf.particles(2,:,:))','y');
h2 = plot(x(2,:));
h3 = plot(pf.xh(2,:),'k');
grid on
title('Y Estimation')
xlabel('Samples')
ylabel('Position (dm)')
legend([h1(1) h2 h3],'Particles','X','X_h')

subplot(3,1,3)
hold on
h1 = plot(squeeze(pf.particles(3,:,:))','y');
h2 = plot(x(3,:));
h3 = plot(pf.xh(3,:),'k');
grid on
title('Theta Estimation')
xlabel('Samples')
ylabel('Orientation (rad)')
legend([h1(1) h2 h3],'Particles','X','X_h')

%% Estimation Error
figure(4)
subplot(3,1,1)
plot(x(1,:)-pf.xh(1,:))
grid on
title('X Estimation Error')
xlabel('Samples')
ylabel('Error (dm)')
subplot(3,1,2)
plot(x(2,:)-pf.xh(2,:))
grid on
title('Y Estimation Error')
xlabel('Samples')
ylabel('Error (dm)')
subplot(3,1,3)
plot(x(3,:)-pf.xh(3,:))
grid on
title('Theta Estimation Error')
xlabel('Samples')
ylabel('Error (rad)')
