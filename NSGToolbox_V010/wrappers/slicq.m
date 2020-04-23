function [c,g,shift,M,Ls,sl_len,tr_area] = ...
    slicq(f,fmin,fmax,bins,sl_len,tr_area,sr,M,min_win,Qvar)
%SLICQ  Sliced constant-Q/variable-Q transform
%   Usage:  [c,g,shift,M,Ls,sl_len,tr_area] 
%                = slicq(f,fmin,fmax,bins,sl_len,tr_area,sr,M,min_win,Qvar)
%                = slicq(f,fmin,fmax,bins,sl_len,tr_area,sr,M,min_win)
%                = slicq(f,fmin,fmax,bins,sl_len,tr_area,sr,M)
%                = slicq(f,fmin,fmax,bins,sl_len,tr_area,sr)
%                = slicq(f,fmin,fmax,bins,sl_len,tr_area)
%                = slicq(f,fmin,fmax,bins,sl_len)
%                = slicq(f,fmin,fmax,bins)
%                = slicq(f,fmin,fmax)
%                = slicq(f,fmin)
%           c = slicq(...)
%
%   Input parameters:
%         f         : Input signal
%         fmin      : Desired minimum frequency (in Hz)
%         fmax      : Desired maximum frequency (in Hz)
%         bins      : Bins per octave (constant or vector (for VQ))
%         sl_len    : Desired slice length (in samples)
%         tr_area   : Transition area length (in samples, <= sl_len/2)
%         sr        : Sampling rate (in Hz)
%         M         : Desired number of time steps per slice, if set to 
%                     0, a channel vector will be computed (*M must be a 
%                     multiple of 4 or will be set to 4*ceil(M/4))
%         min_win   : Minimum window bandwidth (default 16 samples)
%         Qvar      : Factor varying the bandwidth. Qvar=X leads to a 
%                     Q-factor of Q/X
%   Output parameters:
%         c         : Cell array of coefficients
%         g         : Cell array of analysis filters
%         shift     : Vector of frequency shifts of filters
%         M         : Number of time steps per slice (vector or constant)
%         Ls        : Original signal length
%         sl_len    : Slice length
%         tr_area   : Transition area length
%
%   This is a wrapper function for the sliced constant-Q nonstationary 
%   Gabor transform of the signal f. The signal is smoothly sliced into
%   half-overlap segments of length sl_len weighted by a Tukey window 
%   with transition areas of length tr_area and total length of 
%   sl_len/2 + tr_area. 
%
%   Subsequently, a constant-Q nonstationary Gabor transform with essential
%   frequency range fmin to fmax and bins bins per octave will be
%   applied to each segment using NSGCQWIN with modified Blackman-Harris 
%   windows and NSGTF. 
%
%   The additional parameters are an optional fixed number of time steps 
%   M per slice in each frequency channel and a bandwidth variation 
%   factor Qvar. Setting the minimum support min_win of the filters 
%   used helps in preserving shape and localization of low frequency
%   filters, but may lead to a varying Q-factor in that frequency range.
%   
%   See the help of NSGCQWIN for more information on the constant-Q
%   nonstationary Gabor transform.
%
%   See also:  islicq, nsgtf, nsgcqwin, slicing
%
%   References:
%     G. A. Velasco, N. Holighaus, M. Dörfler, and T. Grill. Constructing an
%     invertible constant-Q transform with non-stationary Gabor frames.
%     Proceedings of DAFX11, Paris, 2011.
%     
%     N. Holighaus, M. Dörfler, G. Velasco, and T. Grill. A framework for
%     invertible, real-time constant-q transforms. Audio, Speech, and
%     Language Processing, IEEE Transactions on, 21(4):775-785, April 2013.
%     
%
%   Url: http://nsg.sourceforge.net/doc/wrappers/slicq.php

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

if size(f,1) == 1
    f = f.';
end

Ls = length(f);

if nargin < 10
    Qvar = 1;
    if nargin < 9
        min_win = 16;
        if nargin < 8
            M = 0;
            if nargin < 7
                sr = 1;
                if nargin < 6
                    if nargin < 5
                        sl_len = 16384;
                        if nargin < 4
                            bins = 12;
                            if nargin < 3
                                fmax = .5;
                                if nargin < 2
                                    error('Too few input arguments');
                                end
                            end
                        end
                    end
                    tr_area = round(sl_len/16);
                end
            end
        end
    end
end

if numel(M) > 1
    warning('Number of channels must be a constant or 0, channel ',...
        'vector will be recomputed');
    M(1) = 0;
end

if M(1) == 0
    [g,shift,M] = nsgcqwin(fmin,fmax,bins,sr,sl_len,'min_win',min_win,...
        'Qvar',Qvar,'bwfac',4,'fractional',1,'winfun','modblackharr');
else
    [g,shift] = nsgcqwin(fmin,fmax,bins,sr,sl_len,'min_win',min_win,...
        'Qvar',Qvar,'bwfac',4,'fractional',1,'winfun','modblackharr');
    M = 4*ceil(M/4);
end
N = length(shift); % The number of filters

%% Compute the slices
f_sliced = slicing(f,sl_len,tr_area,Ls);

%% Compute the CQ of each slice
c = nsgtf(f_sliced,g,shift,M);

%% Rearrange the coefficients such that they the slices are centered

if iscell(c) == 0 % Matrix coefficients
    A = [3*M(1)/4+1:M(1),1:3*M(1)/4];
    B = [M(1)/4+1:M(1),1:M(1)/4];
    
    c(:,:,1:2:end) = c(A,:,1:2:end);
    c(:,:,2:2:end) = c(B,:,2:2:end);
else % Cell array coefficients
    for jj = 1:length(g) % Can this be faster somehow?
        A = [3*M(jj)/4+1:M(jj),1:3*M(jj)/4];
        B = [M(jj)/4+1:M(jj),1:M(jj)/4];
        
        c{jj}(:,1:2:end) = c{jj}(A,1:2:end);
        c{jj}(:,2:2:end) = c{jj}(B,2:2:end);
    end
end

