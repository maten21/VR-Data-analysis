function value = RawDataSizeGB(sessionObject)
%RAWDATASIZEGB Get value for RawDataSizeGB
%   Detailed explanation goes here

% Initialize output value with the default value.
value = nan;                 % Please do not edit this line

% Return default value if no input is given (used during config).
if nargin < 1; return; end	% Please do not edit this line


% Insert your code here:

d = dir(fullfile(sessionObject.DataLocation(1).RootPath, sessionObject.DataLocation(1).Subfolders,  '**\*.*'));

b = 0;
for i = 1:length(d)
    b = b + d(i).bytes;
end

value = b / (1024*1024*1024);


end

