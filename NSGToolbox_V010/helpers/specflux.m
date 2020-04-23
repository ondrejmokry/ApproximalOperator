function [SF,V0] = specflux(f,win_length,tgap)
%SPECFLUX  Spectral flux onset detection function
%   Usage: [SF,V0] = specflux(f,win_length,tgap)
%          SF = specflux(f,win_length,tgap)
%
%   Input parameters: 
%         f         : Input signal
%         win_length: Desired window length for the STFT
%         tgap      : Time step for the STFT
%   Output parameters:
%         SF        : Spectral flux of f*
%         V0        : STFT coefficients of f*
% 
%   This is a helper function for ONSETDET and not meant to
%   be used individually.
%
%   Computes the spectral flux onset-detection function
%   of f with a Hann window of length win_length. 
%   The STFT is taken with time shift parameter tgap*
%   and win_length frequency channels.
%
%   Externals: COMP_DGT_FB (LTFAT routine, included in NSGToolbox V0.1.0 
%              and higher)
%
%   See also:  onsetdet
%
%   References:
%     S. Dixon. Onset detection revisited. In Proceedings of the 9th
%     International Conference on Digital Audio Effects, volume 120, pages
%     133-137, 2006.
%     
%
%   Url: http://nsg.sourceforge.net/doc/helpers/specflux.php

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
% Date: 26.04.13

% Check input arguments

if nargin < 3
    error('Not enough input arguments');
end

% Compute the Gabor transform (sampled STFT) of f

win=winfuns('hann',win_length);

if size(f,2) > 1
    f = f.';
end

V0 = comp_dgt_fb(f,win,tgap,win_length);

% Compute the spectral flux 

VV = abs(V0);
VV = max(VV-circshift(VV,[0,1]),0);

SF = sum(VV);

% Normalize

SF = SF-mean(SF);
SF = SF./std(SF);
