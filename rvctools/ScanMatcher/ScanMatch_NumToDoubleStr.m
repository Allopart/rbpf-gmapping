function seq = ScanMatch_NumToDoubleStr(num, modulus)
%SCANMATCH_NUMTODOUBLESTR transform an array of numbers into a sequence of
% double characters (each number is represented by a residue of two letter) in
% function of the modulus. The modulus is the length of the
% alphabet used for the second letter of a residue. 
% i.e.  seq = ScanMatch_NumToDoubleStr([1 2 27], 26) % modulus = 26 (all the alphabet is used)
%  returns aAaBBa
% i.e.  seq = ScanMatch_NumToDoubleStr([1 2 27], 10) % modulus = 26 (all the alphabet is used)
%  returns aAaBcG
%
% The modulus can be any number between 1 and 26.
% The Input numbers have to be between 1 and 676 (which is the maximum you 
% can code with two characters: 676 is xX with a modulus of 26 )
%
%   Part of the ScanMatch toolbox
%   Written by Filipe Cristino 
%   $Version: 1.00 $  $Date: 10/09/2009

% ---- Check Modulus ----
if modulus > 26 || modulus < 1
    error('ScanMatch:WrongInput', 'The modulus has to be between 1 and 26!')
end

% ---- Check input numbers ----
max_num = 26 * modulus; % Max number you can convert in function of the modulus
if sum(num <1 | num > max_num)
    error('ScanMatch:WrongInput', 'The input number can not be less than 1 or more than (26 * modulus)')
end

% ---- substract one to the array as the conversion works from zero
num = num - 1; 

% ---- Start conversion ----
num_l = length(num);
ind= 1;
for i=1:num_l
    seq(ind) = lower(char((fix(num(i)/modulus))+65));
    seq(ind+1) = char(rem(num(i),modulus)+65);
    ind = ind+2;
end