function MF_startup()
% MF_STARTUP  Startup file for MultiFrontal
%   MAKE adds paths of the MultiFrontal to Matlab.

%   Copyright (c) 2015 Yingzhou Li, Stanford University

file_path = mfilename('fullpath');
tmp = strfind(file_path,'MF_startup');
file_path = file_path(1:(tmp(end)-1));

addpath(genpath([file_path 'src']));

end
