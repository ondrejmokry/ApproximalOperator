function fr = invwvlttrans(c,g,shift,M,tgtfl,fb,Ls)
%INVWVLTTRANS  Wavelet frame synthesis
%   Usage:  fr = invwvlttrans(c,g,shift,M,tgtfl,fb,Ls)
%           fr = invwvlttrans(c,g,shift,M,tgtfl,fb)
%           fr = invwvlttrans(c,g,shift,M,tgtfl,Ls)
%           fr = invwvlttrans(c,g,shift,M,tgtfl)
%           fr = invwvlttrans(c,g,shift,M,Ls)
%           fr = invwvlttrans(c,g,shift,M)
% 
%   Input parameters: 
%         c         : Cell array of Wavelet coefficients
%         g         : Cell array of Fourier transforms of the analysis 
%                     Wavelets
%         shift     : Vector of frequency shifts
%         M         : Number of time steps
%         Ls        : Original signal length
%         fb	    : Frame bounds (vector)
%         tgtfl     : Tightflag (1 if frame is tight)
%   Output parameters:
%         fr        : Reconstructed signal
%
%   This is a wrapper function for the inverse painless Wavelet transform
%   via nonstationary Gabor filterbank. Given the cell array c and the 
%   painless Wavelet frame g, shift, M, the corresponding dual frame 
%   is computed and the corresponding inverse Wavelet transform is 
%   performed.
%
%   If the original signal length Ls is specified, the synthesized signal
%   will be truncated to length Ls. If the parameters tgtfl and fb*
%   are given, the system g, shift, M is assumed to be a tight frame
%   and synthesis is performed using the original system.
%
%   See also:  wvlttrans, nsigtf, nsdual
%
%   References:
%     P. Balazs, M. Dörfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
%     Theory, implementation and applications of nonstationary Gabor Frames.
%     J. Comput. Appl. Math., 236(6):1481-1496, 2011.
%     
%
%   Url: http://nsg.sourceforge.net/doc/wrappers/invwvlttrans.php

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

if nargin < 7
    if numel(fb) == 1 && tgtfl
        Ls = fb;
        fb = [1 1]
    elseif numel(fb) == 1
        Ls = fb;
        tgtfl = 0;
    end
    if nargin < 6
        if nargin < 4
            error('Too few input arguments');
        end
        if exist('tgtfl','var')
            if tgtfl > 1
                Ls = tgtfl;
                tgtfl = 0;
            end
        else
            tgtfl = 0;
        end
    end
end

% Compute the dual frame
if tgtfl
    gd =g;
else
    gd = nsdual(g,shift,M);
end

% Inverse Wavelet transform
if nargin == 4
    fr = nsigtf(c,gd,shift);
else
    fr = nsigtf(c,gd,shift,Ls);
end

if tgtfl
    fr = fr/fb(1);
end
