function Ok = ScanMatch_CheckStructure(ScanMatchInfo)
% SCANMATCH_CHECKSTRUCTURE checks if the ScanMatchInfo structure contains all the
% required fields. If the structure is complete it returns 1 if not an
% error message is return and the execution of the program is stopped.
%
%   Part of the ScanMatch toolbox
%   Written by Filipe Cristino 
%   $Version: 1.00 $  $Date: 10/09/2009


OkFields = isfield(ScanMatchInfo, { 'Xres','Yres','Xbin', 'Ybin',...
        'RoiModulus','Threshold','mask','SubMatrix', 'GapValue', 'TempBin'});
    
if sum(OkFields) > 9
    Ok = 1;
else
    error('ScanMatch:MissingFields', 'The ScanMatchInfo is missing %d fields', 10 - sum(OkFields))
end