function [ lines ] = RansacLines( laserScanCartesian, RNSC)
% [lines] = RANSACLINES(laserScanCartesian, parameters)
% This function extracts a set of lines from the given set of points in the
% cartesian coordinates. The algorithm used is random sample consensus
% (RANSAC).

noOfPoints = size(laserScanCartesian,2);

%% Parameters
maxNoOfLines = RNSC.MaxIter; % The maximum number of line extraction iterations
noOfRandomCouples = RNSC.Couples; % The number of random point couples to try before choosing the best candidate
distThreshold= RNSC.Threshold; % The distance threshold determining whether a point is supporting a line
minLineSupport= RNSC.MinLineSupport; % The minimum number of points to support a line for the line to be accepted
minNoOfPoints= RNSC.MinNoOfPoints; % The minimum number of points to continue line extraction iterations

ransacLinesDemo = 0; % setting this to 1 activates some plotting for demonstration

%% Pre Processing
lines = zeros(2,maxNoOfLines);
linesInd = 1;
%% Main Loop
for i=1:maxNoOfLines
    if(noOfPoints<minNoOfPoints) %Check if enough points are left
        break
    end
    
    maxSupport = 0;
    for j=1:noOfRandomCouples % Find a candidate line by trying a number of random point couples
        p1 = laserScanCartesian(:,randint(1,1,[1,noOfPoints]));
        p2 = laserScanCartesian(:,randint(1,1,[1,noOfPoints]));
        while(all(p1==p2))
            p2 = laserScanCartesian(:,randint(1,1,[1,noOfPoints]));
        end
        points=[p1,p2];
        candLine = lsqLine(points);
        dists=cos(candLine(1))*laserScanCartesian(1,1:noOfPoints)+sin(candLine(1))*laserScanCartesian(2,1:noOfPoints)-candLine(2);
        oldAdmit = abs(dists)<distThreshold;
        lineSupport = sum(oldAdmit);
        if(lineSupport>maxSupport)
            maxSupport=lineSupport;
            bestAdmit = oldAdmit;
            bestLine = candLine;
        end
    end
    candLine = bestLine; %CAREFUL I WAS COMMENTED OUT!
    oldAdmit = bestAdmit;
    if(lineSupport>minLineSupport) %If the resulting line is acceptable, iterate on it to fit better
        while(1)
            points = laserScanCartesian(:,oldAdmit);
            candLine = lsqLine(points);
            dists=cos(candLine(1))*laserScanCartesian(1,1:noOfPoints)+sin(candLine(1))*laserScanCartesian(2,1:noOfPoints)-candLine(2);
            admit = abs(dists)<distThreshold;
            if(all(oldAdmit==admit))
                break;
            end
            oldAdmit=admit;
        end
    lines(:,linesInd)=candLine;
    if(ransacLinesDemo)
        plot(laserScanCartesian(1,:)',laserScanCartesian(2,:),'b')
        hold on
        plot(laserScanCartesian(1,admit)',laserScanCartesian(2,admit),'.r')
        hold off
        pause
    end
        linesInd=linesInd+1;
        laserScanCartesian = laserScanCartesian(:,~admit);
        noOfPoints = size(laserScanCartesian,2);
    end
end

lines = lines(:,1:linesInd-1);
lines = lines(:,~any(isnan(lines),1));

end