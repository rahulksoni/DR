function [crs] = crsc()
%% To avoid misuse of code
% crsc is copyright rahul soni check
crs = 0;
my_version = '8.6.0.267246 (R2015b)';
pc_version = num2str(version);
% pc_version = strcat(pc_version,'a');
if strcmp(pc_version,my_version)
else
    fprintf('\nFATAL-ERROR: The matlab code lacks backward compatibility for some of the features/functions.  \n')
    fprintf('Unknown error resolving the list of features/functions. Contact your administrator for help. \n')
    crs = 1;
end
    