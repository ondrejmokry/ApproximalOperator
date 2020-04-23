function [g,shift,M] = nsgsclwin(positions,short,max_win,Ls)
%NSGSCLWIN  Scale-frame dictionary generator
%   Usage:  [g,shift,M] = nsgsclwin(positions,short,max_win,Ls)
%
%   Input parameters:
%         positions : A vector specifying time positions around which
%                     shortest windows are desired
%         short     : Shortest admissible window length (in samples)
%         max_win   : Maximum number of different window sizes
%         Ls        : Length of the signal to be analyzed (in samples)
%
%   Output parameters:
%         g         : Cell array of window functions
%         shift     : Vector of time shifts
%         M         : Vector of window lengths
%
%   Create nonstationary Gabor frame from a sequence of positions, e.g. 
%   an onset sequence. This routine builds scale-type frames with fixed 
%   scale parameter Q = 2 and overlap parameters O_1 = 2/3, 
%   O_2 = 1/3. Currently, Hann windows of varying length are used as 
%   prototype. 
%
%   The routine places short windows on the points specified by the
%   sequence, while the space between 2 points is spanned by a sequence of
%   that smoothly expand to longer lengths and then shrink back to a short
%   length. We call those short sequences short->long->short building 
%   blocks.
%
%   The final window sequence is obtained by concatenating the sequence of
%   building blocks. 
% 
%   For more details, please check the reference.
%
%   Notes: 
%
%   For the parameters Q, O_1 and O_2 as specified above being valid,
%   the input parameter short must be a multiple of 6, or non integer 
%   shifts might occur.
% 
%   The first value of positions should always be 1 to cover the complete
%   time axis. A valid position sequence can be obtained, e.g. from
%   ONSETDET.
%
%   See also:  nsgt, nsgt_real, onsetnsgt, invonsetnsgt, onsetdet
%
%   References:
%     P. Balazs, M. Dörfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
%     Theory, implementation and applications of nonstationary Gabor Frames.
%     J. Comput. Appl. Math., 236(6):1481-1496, 2011.
%     
%
%   Url: http://nsg.sourceforge.net/doc/generators/nsgsclwin.php

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
% Date: 23.04.13

if nargin < 4
    error('Not enough input arguments');
end
if rem(short,6) ~= 0
    error('Shortest admissible window length must be a multiple of 6');
end

ss_3 = short/3;

positions = [positions;Ls]; % The last window should be 'shortest', so
% the length of the signal is added to the
% onset vector.

% 'blocks' specifies the step size between positions.

blocks = positions(2:end)-positions(1:end-1) - 2*ss_3;

if numel(find(blocks < 2*ss_3)) > 0
    error('Positions vector and shortest admissible window size ',...
        'incompatible');
end

NN = length(blocks);

% 'dists(n)' holds the minimal length of 'short->long->short'
% building blocks with 'n' different windows

dists = (7*2.^(0:max_win-1) - 5)*ss_3;

% 'm0(n)' will be the number of different windows used in the
% n-th 'short->long->short' building block

%m0 =  max_win*ones(NN,1);

% 'A(n,m)' will be the number of fill-in windows of scale
% 'm' needed to stretch the n-th building block to the
% desired size

A = zeros(NN,max_win);

% Determine, from 'blocks', the correct values of 'm0'
% and 'A' for the construction of the frame

m0 = arrayfun(@(x) find( dists <= x,1,'last'),blocks);
blocks = arrayfun(@(x,y) x-dists(y),blocks,m0);

idx = sub2ind(size(A),(1:NN)',m0);
A(idx) = floor(blocks./(ss_3*2.^m0));

blocks = arrayfun(@rem,blocks,ss_3*2.^m0);

fillin = mod(cumsum(blocks)-1,2*ss_3)+1;

blocks(1:end) = blocks(1:end)-fillin(1:end)+[0;fillin(1:end-1)];
blocks(1) = blocks(1)+2*ss_3;

temp = cell2mat(arrayfun(@(x) de2bi(x/ss_3,max_win+1),...
    blocks,'UniformOutput',0));

A = A + temp(:,2:end);

shift = zeros(NN,1); % Initialize shift vector
win_size = shift;
n = 1;

% Determine, from 'm0' and 'A', the window lengths
% and time shifts corresponding to time step 'n'.

for k = 1:NN
    nk = n+A(k,1)+1;
    win_size(n:nk-1) = short;
    shift(n:nk-1) = 2*ss_3;
    n = nk;
    for l = 2:m0(k)-1
        nk = n+A(k,l)+1;
        Q = 2^(l-1);
        win_size(n:nk-1) = short*Q;
        shift(n) = 5/4*ss_3*Q;
        shift(n+1:nk-1) = 2*ss_3*Q;
        n = nk;
    end
    l = m0(k);
    Q = 2^(l-1);
    if l > 1
        win_size(n) = short*Q;
        shift(n) = 5/4*ss_3*Q;
        n = n+1;
    end
    nk = n+A(k,l)+1;
    win_size(n:nk-1) = short*Q;
    shift(n:nk-1) = 2*ss_3*Q;
    n = nk;
    nk = n+m0(k)-1;
    Q = 2.^(m0(k)-2:-1:0)';
    win_size(n:nk-1) = short*Q;
    shift(n:nk-1) = 5/2*ss_3*Q;
    n = nk;
end

% Calculate the windows corresponding to the time
% steps (can be improved in various ways)

M = win_size;

g = arrayfun(@(x,y) winfuns('hann',x)/sqrt(y),win_size,M,'UniformOutput',0);
