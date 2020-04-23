function nsgt_startup()
%NSGT_STARTUP  Set paths for using NSGToolbox
%   Usage: nsgt_startup()
%
%   This script file adds NSGToolbox folders to the MATLAB path.
%
%   Url: http://nsg.sourceforge.net/doc/nsgt_startup.php

% Copyright (C) 2013 Nicki Holighaus.
% This file is part of NSGToolbox version 0.1.0
% 
% This work is licensed under the Creative Commons 
% Attribution-NonCommercial-ShareAlike 3.0 Unported 
% License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-sa/3.0/ 
% or send a letter to 
% Creative Commons, 444 Castro Street, Suite 900, 
% Mountain View, California, 94041, USA.

% Author: Nicki Holighaus
% Date: 20.03.13

basepath = which('nsgt_startup');
basepath=basepath(1:end-15);

addpath(basepath);

d = dir(basepath);

for ii=1:length(d)
  
  % We only look for directories
  if ~d(ii).isdir
    continue;
  end;
  
  % Skip the default directories . and ..
  if (d(ii).name(1)=='.')
    continue;
  end;
  
  addpath([basepath,filesep,d(ii).name]);    
end

basepath = [basepath,filesep,'wrappers'];
d = dir(basepath);

for ii=1:length(d)
  
  % We only look for directories
  if ~d(ii).isdir
    continue;
  end;
  
  % Skip the default directories . and ..
  if (d(ii).name(1)=='.')
    continue;
  end;
  
  addpath([basepath,filesep,d(ii).name]);    
end
