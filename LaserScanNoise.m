
% -------------------------------------------------------------------------
% laserscan.m
%
% Filen returnerer de synlige punkter for en laserscanner med en max. 
% afstand på 15m, set fra en given position og retning i et kendt område. 
% Resultatet returneres i form af et 2x361-array med henholdsvis vinkler og
% afstande.
%
% Skrevet af Christian Overvad, s031914, og Kasper Strange, s031924.
% 
% Sidst ændret 06-06-2006
%  
%
% Returns the visible points for a laser scanner with a max distance at 
% maxDistance meters, seen from a given position and direction in a 
% well-known area being (x,y,theta) the CURRENT GLOBAL POSITION OF THE
% LASER SCANNER. Scan width is restricted to 180 degrees.
% The result is returned in an 2xN-array with respectively angles and 
% distances.
% 
% Edited by: Rafael Olmos Galve s071150 at DTU, Denmark - rafa_olmos@hotmail.com
% Changes done: -Max laser scanner measurement distance as a parameter.
%               -Resolution is not restricted anymore, it is specified as
%               parameter.
% Last change 12/05/2008.
%
% Rewritten by Nils Axel Andersen 23-2-2011
% -------------------------------------------------------------------------

function scan = LaserScanNoise(pose, lines, LRS, sigma)
%Extract Struct
maxDistance     = LRS.MaxDistance;
resol           = LRS.Resolution;
field_of_view   = LRS.FoV;

%Extract Pose
x = pose(1);
y = pose(2);
theta = pose(3);


% Number of scanning lines (deduced from scan width and resolution)
num_of_scanning_lines=floor(field_of_view/resol)+1;
resolrad=resol*pi/180.0;
% -------------------------------------------------------------------------
% Function [M,N] = SIZE(X) for matrix X, returns the number of rows(M) and 
% columns(N) in X as separate output variables.
% -------------------------------------------------------------------------

[~, no_of_lines] = size(lines); % Totalt antal linier - Total num. of lines

% ---------------------------------------------------
% Præ-allokering af arrays - Preallocation of arrays
% ---------------------------------------------------

% Global system coordinates
x_start(1:no_of_lines) = 0;
y_start(1:no_of_lines) = 0;
x_end(1:no_of_lines) = 0;
y_end(1:no_of_lines) = 0;

% Laser scanner local system coordinates
trans_x_start(1:no_of_lines) = 0;
trans_y_start(1:no_of_lines) = 0;
trans_x_end(1:no_of_lines) = 0;
trans_y_end(1:no_of_lines) = 0;

b(1:no_of_lines) = 0;
a(1:no_of_lines) = 0;
scan(1:2,1:no_of_lines) = 0;

% --------------------------------------------------
% Start- og slutværdier for liniernes endepunkterne.
% Start and end values for the lines ending points.
% --------------------------------------------------

for i = 1:no_of_lines
    x_start(i) = lines(1,i);
    y_start(i) = lines(2,i);
    x_end(i)   = lines(3,i);
    y_end(i)   = lines(4,i);
    
% ---------------------------------------------------------------------
% Transformation af linierne til laserscannerens koordinatsystem.
% Transformation of the lines to the laser scanner's coordinate system.
% ---------------------------------------------------------------------   

    % Linierne konverteres til nye akser.
    % Lines are converted to the new coordinate system.
    trans_x_start(i)=(x_start(i)-x)*cos(theta) + (y_start(i)-y)*sin(theta);
    trans_y_start(i)=(y_start(i)-y)*cos(theta) - (x_start(i)-x)*sin(theta); 
    trans_x_end(i) = (x_end(i)-x)*cos(theta) + (y_end(i)-y)*sin(theta);
    trans_y_end(i) = (y_end(i)-y)*cos(theta) - (x_end(i)-x)*sin(theta);
    
    % Den mindste x-værdi sættes til x_start(i).
    % Starting and ending points are swapped if the x value of the 
    % starting point is bigger than the x value of the ending point.
    if trans_x_start(i) > trans_x_end(i)
        trans_x_temp = trans_x_start(i);
        trans_x_start(i) = trans_x_end(i);
        trans_x_end(i) = trans_x_temp;
        trans_y_temp = trans_y_start(i);
        trans_y_start(i) = trans_y_end(i);
        trans_y_end(i) = trans_y_temp;
    end;

% -----------------------------------------------------------------------
% The line equations are calculated 
%  a*x+ b*y=c
% -----------------------------------------------------------------------
  a(i)=trans_y_start(i)-trans_y_end(i);
  b(i)=trans_x_end(i)-trans_x_start(i);
  l=sqrt(a(i)^2+b(i)^2);
  a(i)=a(i)/l;
  b(i)=b(i)/l;
  c(i)=a(i)*trans_x_start(i)+b(i)*trans_y_start(i);
end;

% -------------------------------------------------------------------
% Konvertering af det, som laserscanneren ser til polære koordinater.
% Conversion of what the laser scanner sees to polar coordinates.
% -------------------------------------------------------------------
% For each laser scanner angle
for i = 1:num_of_scanning_lines
    % Laserscanerens maksimale måleafstand.
    % Laser scanner maximum measured distance maxDistance(in meters)
    % Closest distance from the laser scanner.
    min_dist = maxDistance;
    % Aktuelle vinkel laserscanner.
    % Current laser scanner angle 
    phi =(i-1)*resolrad-field_of_view/2.0*pi/180.0;   
    % linierne gennemgås for at finde deres afstand til laserscanneren
    % Find distance from the lines to laser scanner in the current angle
    cosphi=cos(phi);
    sinphi=sin(phi);
    for j = 1:no_of_lines
     temp=a(j)*cosphi+b(j)*sinphi;
     if ( abs(temp)>1e-6)
        t=c(j)/temp;
        if (t>0 && t<min_dist)
           if (abs(trans_x_start(j)-trans_x_end(j))>1e-6 )
               if(t*cosphi < trans_x_end(j) && t*cosphi > trans_x_start(j))
                 min_dist=t;
               end
           else
             if (trans_y_end(j) > trans_y_start(j) )
                if(t*sinphi < trans_y_end(j) && t*sinphi > trans_y_start(j))
                  min_dist=t; 
                end
             else
                if(t*sinphi > trans_y_end(j) && t*sinphi < trans_y_start(j))
                  min_dist=t; 
                end
             end
        end
     end
     end    
    end
%     if min_dist == maxDistance
%         %scan(:,i) = nan;
%         scan(1,i) = phi;
%         scan(2,i) = min_dist;   
%     else
    % The polar coordinates returned
    scan(1,i) = phi + sigma(1)*randn;
    scan(2,i) = min_dist + sigma(2)*randn;
    %end
end
    %scan(:,isnan(scan(1,:)))=[];