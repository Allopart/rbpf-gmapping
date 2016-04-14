function res = ScanMatch_DoubleStrToNum(str, modulus)
% SCANMATCH_DOUBLESTRTONUM transforms a string sequence of double characters 
% into a single number in function of the modulus used. The modulus is the 
% lenght of the alphabet used for the second letter of a residue. 
%  i.e. res = ScanMatch_DoubleStrToNum('aAaBbA', 26)% Modulus = 26 (all the alphabet is used)
%   returns     1 2 27
%  or   res = ScanMatch_DoubleStrToNum('aAaBbA', 10)% Modulus = 10 (only aA to aJ is used)
%   returns     1 2 11
%
% The modulus can be any number between 1 and 26.
% The Input numbers have to be between 1 and 676 (which is the maximum you can code with two:  
%    676 is xX with a modulus of 26 )
%
%   Part of the ScanMatch toolbox
%   Written by Filipe Cristino 
%   $Version: 1.00 $  $Date: 10/09/2009

% ---- Check Modulus ----
if modulus > 26 || modulus < 1
    error('ScanMatch:WrongInput', 'The modulus has to be between 1 and 26!')
end

% Check needed for the second letter of a residue in function of teh
% modulus

% Upper case the string
str = upper(str);

% lengh of string
str_l = length(str);

% check if string is even
if mod(str_l,2) ~= 0
    error('ScanMatch:WrongInput','The number of characters in the input string need to be even...')
end

% ---- do the processing two by two letters ----
ind=1;
res = zeros(1,str_l/2);
for i=1:2:str_l
    num_str = double(str(i:i+1)) - 64;
    res(1,ind) = (num_str(1)-1)* modulus + num_str(2) - 1;
    ind = ind+1;
end

% Add one to the final array as the conversion works from zero
res = res + 1;

end


