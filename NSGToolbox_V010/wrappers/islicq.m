function [fr,gd] = islicq(c,g,shift,M,Ls,sl_len,tr_area)
%ISLICQ  Sliced constant-Q/variable-Q synthesis
%   Usage:  [fr,gd] = islicq(c,g,shift,M,Ls,sl_len,tr_area)
%           fr = islicq(c,g,shift,M,Ls,sl_len,tr_area)
%
%   Input parameters:
%         c         : Cell array of coefficients
%         g         : Cell array of analysis filters
%         shift     : Vector of frequency shifts of filters
%         M         : Number of channels (vector/scalar)
%         Ls        : Original signal length
%         sl_len    : Slice length
%         tr_area   : Transition area length
%   Output parameters:
%         fr        : Reconstructed signal
%         gd        : Cell array of synthesis windows
%
%   This is a wrapper function for the inverse sliced constant-Q 
%   nonstationary Gabor transform. Given an array of coefficients c 
%   corresponding to a sliced nonstationary Gabor system g, shift, M 
%   with slice length sl_len, this function computes the corresponding
%   synthesis operation. 
%
%   That is, NSDUAL is used to compute the canonical dual frame of the
%   system g, shift, M if possible. Afterwards, a sliced signal is
%   synthesized using this dual system and the function NSIGTF. Finally
%   the synthesized signal is unsliced by the helper function UNSLICING.
%   
%   If the output of SLICQ is used as input for this function, and the
%   system g, shift, M used is a painless frame, then the originally
%   analyzed function is reconstructed perfectly.
%
%   See the help of NSGCQWIN for more information on the constant-Q
%   nonstationary Gabor transform.
%
%   See also:  slicq, nsigtf, nsdual, unslicing
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
%   Url: http://nsg.sourceforge.net/doc/wrappers/islicq.php

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
    error('Not enough input arguments');
end

% Before reconstruction, the coefficients have to be rearranged

if iscell(c) == 0 % Matrix coefficients
    slices = size(c,3);
    A = [3*M(1)/4+1:M(1),1:3*M(1)/4];
    B = [M(1)/4+1:M(1),1:M(1)/4];
    
    c(:,:,1:2:slices) = c(B,:,1:2:slices);
    c(:,:,2:2:slices) = c(A,:,2:2:slices);
else % Cell array coefficients
    slices = size(c{1},2);
    for jj = 1:length(g) % Can this be faster somehow?
        A = [3*M(jj)/4+1:M(jj),1:3*M(jj)/4];
        B = [M(jj)/4+1:M(jj),1:M(jj)/4];
        
        c{jj}(:,1:2:end) = c{jj}(B,1:2:end);
        c{jj}(:,2:2:end) = c{jj}(A,2:2:end);
    end
end

gd = nsdual(g,shift,M); % Compute the dual frame

% The actual transform
fr = isliCQ_int(c,gd,shift,Ls,sl_len,tr_area,slices);

end

%% Internal computation of the transform

function fr = isliCQ_int(c,gd,shift,Ls,sl_len,tr_area,slices)

%% Compute the inverse CQ of each slice
frec_sliced = nsigtf(c,gd,shift,sl_len);
%% Glue the parts back together
fr = unslicing(real(frec_sliced),sl_len,tr_area,slices);
fr = fr(1:Ls);

end

