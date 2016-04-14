function [score path step F] = ScanMatch_nwAlgo(intseq1,m,intseq2,n,ScoringMatrix,gap)
%SCANMATCH_NWALGO Standard Needleman-Wunsch algorithm implementation
% This is originally from the Matlab Bioinformatics toolbox but was partially 
% changed to work with ScanMatch  

%   Part of the ScanMatch toolbox
%   From the BioInformatics Toolbox, Modified by Filipe Cristino 
%   $Version: 1.00 $  $Date: 10/09/2009

% set up storage for dynamic programming matrix
F = zeros(n+1,m+1);
F(2:end,1) = gap * (1:n)';
F(1,2:end) = gap * (1:m);

% and for the back tracing matrix
pointer= repmat(uint8(4),n+1,m+1);
pointer(:,1) = 2;  % up
pointer(1,1) = 1;  


% initialize buffers to the first column
ptr = pointer(:,2); % ptr(1) is always 4
currentFColumn = F(:,1);

% main loop runs through the matrix looking for maximal scores
for outer = 2:m+1

    % score current column
    scoredMatchColumn = ScoringMatrix(intseq2,intseq1(outer-1));
    % grab the data from the matrices and initialize some values
    lastFColumn    = currentFColumn;
    currentFColumn = F(:,outer);
    best = currentFColumn(1);
    
    for inner = 2:n+1
        % score the three options
        up       = best + gap;
        left     = lastFColumn(inner) + gap;
        diagonal = lastFColumn(inner-1) + scoredMatchColumn(inner-1);

        % max could be used here but it is quicker to use if statements
        if up > left
            best = up;
            pos = 2;
        else
            best = left;
            pos = 4;
        end

        if diagonal >= best
            best = diagonal;
            ptr(inner) = 1;
        else
            ptr(inner) = pos;
        end
        currentFColumn(inner) = best;

    end % inner
    % put back updated columns
    F(:,outer)   = currentFColumn;
    % save columns of pointers
    pointer(:,outer)  = ptr;
end % outer

% Find the best route throught the scoring matrix
i = n+1; j = m+1;
path = zeros(n+m,2);
step = 1;

while (i > 1 || j > 1)
    
    switch pointer(i,j)
        case 1 % diagonal only
            j = j - 1;
            i = i - 1;
            path(step,:) = [j,i];
        case 2 % up only
            i = i - 1;
            path(step,2) = i;
        case 4 % left only
            j = j - 1;
            path(step,1) = j;
        case 6 % up or left --> up (favors gaps in seq2)
            j = j - 1;
            path(step,1) = j;
        otherwise %3 diagonal or up   --> diagonal (favors no gaps)
            %4 diagonal or left       --> diagonal (favors no gaps)
            %7 diagonal or left or up --> diagonal (favors no gaps)
            j = j - 1;
            i = i - 1;
            path(step,:) = [j,i];
    end
    step = step +1;
end

score =  max(F(n+1,m+1,:));

    
