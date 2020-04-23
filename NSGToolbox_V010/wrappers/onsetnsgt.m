function [c,g,shift,M,Ls] = onsetnsgt(f,thre,short,max_win,win_length)
%ONSETNSGT  Onset-based nonstationary Gabor transform
%   Usage:  [c,g,shift,M,Ls] = onsetnsgt(f,thre,short,max_win,win_length)
%           [c,g,shift,M,Ls] = onsetnsgt(f,thre,short,max_win)
%           [c,g,shift,M,Ls] = onsetnsgt(f,thre,short)
%           [c,g,shift,M,Ls] = onsetnsgt(f,thre)
%           c = onsetnsgt(...)
%
%   Input parameters: 
%         f         : The signal to be analyzed (single channel only)
%         thre      : Peak-picking threshold
%         short     : Shortest allowed window length
%         max_win   : Maximum number of different windows used
%         win_length: Window length for the onset STFT analysis 
%   Output parameters: 
%         c         : Cell array of transform coefficients
%         g         : Cell array of analysis windows
%         shift     : Vector of time shifts
%         M         : Number of frequency channels
%         Ls        : Original signal length
%
%   This is a wrapper function for the scaleframe nonstationary Gabor
%   transforms with onset detection based adaptation. Given a signal f,
%   this wrapper computes the spectral flux onset detection function based
%   on a regular discrete Gabor transform with redundancy 16 using a Hann
%   window of length win_length. A simple peakpicking algorithm
%   determines the significant maxima in the spectral flux function. Those
%   are assumed to be the onsets in f.
%
%   From this onset sequence, a scaleframe nonstationary Gabor system will
%   be constructed and the corresponding analysis performed by NSGT_REAL.
%  
%   Note: The current wrapper only supports the threshold parameter thre 
%   of the onset detection algorithm. To obtain optimal results, the 
%   remaining parameters need to be fine tuned as well. An experienced user
%   should use ONSETDET, NSGSCLWIN and NSGT_REAL on separately
%   instead. Also see the Onset How-To included in the toolbox.
%
%   See also:  invonsetnsgt, nsgt_real, nsgsclwin, onsetdet
%
%   References:
%     S. Dixon. Onset detection revisited. In Proceedings of the 9th
%     International Conference on Digital Audio Effects, volume 120, pages
%     133-137, 2006.
%     
%     P. Balazs, M. Dörfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
%     Theory, implementation and applications of nonstationary Gabor Frames.
%     J. Comput. Appl. Math., 236(6):1481-1496, 2011.
%     
%
%   Url: http://nsg.sourceforge.net/doc/wrappers/onsetnsgt.php

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
% Date: 25.04.13

if nargin < 5
    win_length = 2048;
    if nargin < 4
        max_win = 10;
        if nargin < 3
            short = 192;
            if nargin < 2
                error('%s: Not enough input arguments',upper(mfilename));
            end
        end
    end
end

if min(size(f)) > 1
    error('%s: Multichannel signals are not supported',upper(mfilename));
end

Ls = length(f);

positions = onsetdet(f,win_length,thre);

blocks = diff(positions);
idx = find(blocks < 4/3*short);
while numel(idx) > 0
    positions(1+idx) = [];
    blocks = diff(positions);
    idx = find(blocks < 4/3*short);
end

[g,shift,M] = nsgsclwin(positions,short,max_win,Ls);

c = nsgt_real(f,g,shift,M);
