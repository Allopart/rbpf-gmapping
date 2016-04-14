function seq_str = ScanMatch_FixationToSequence(data, ScanMatchInfo)
% SCANMATCH_FIXATIONTOSEQUENCE transforms an array of fixation locations and 
% times into double character string sequences.. The input ScanMatchInfo should
% hold the necessary information to perform this. If ScanMatchInfo.TempBin 
% is set to zero then no temporal binning is performed (a single residue 
% per fixation will be used)  
%
%   Part of the ScanMatch toolbox
%   Written by Filipe Cristino 
%   $Version: 1.00 $  $Date: 10/09/2009


% Check if the ScanMatchInfo structure is Ok
Ok = ScanMatch_CheckStructure(ScanMatchInfo);

% Check if the data is at the right format
dataS = size(data);
if dataS(1) == 0
    error('ScanMatch:NoData', 'ScanMatch_FixationToSequence: No Input data!')
end

if dataS(2) > 3
    error('ScanMatch:WrongData', 'ScanMatch_FixationToSequence: Wrong size input data 3 columns max!')
end
    
% ---- Clean the data ----
% Any negative value will be set to 1
data = uint16(data);
data(data == 0)=1;

%Fixations outside the screen resolution will be set to the screen resolution
data = double(data);
data(data(:,1)>ScanMatchInfo.Xres,1) = ScanMatchInfo.Xres-1;
data(data(:,2)>ScanMatchInfo.Yres,2) = ScanMatchInfo.Yres-1;

% ---- Get eye movement sequences ----
subs = sub2ind(size(ScanMatchInfo.mask),data(:,2),data(:,1));
seq_num = ScanMatchInfo.mask(subs)';

% ---- Temporal binning if needed ----
if ScanMatchInfo.TempBin ~= 0 % do this if TempBin is set
    % Check if fixation times are available
    S = size(data);
    if S(2)==3 % fixation times available
        fix_time = round(data(:,3) / ScanMatchInfo.TempBin);
        str = [];
        for i=2:size(data,1)
            str = [str repmat(seq_num(i),1,fix_time(i))];
        end
        seq_num = str;
    else % throw error
        error('ScanMatch:NoFixationtimes', 'No Fixation times found. Set ScanMatchInfo.TempBin to 0 or add fixation times to data')
    end
end

% ---- Create the string sequences ----
 seq_str = ScanMatch_NumToDoubleStr(seq_num,ScanMatchInfo.RoiModulus);

    
end

