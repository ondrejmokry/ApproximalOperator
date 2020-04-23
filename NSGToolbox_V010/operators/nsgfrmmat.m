function S = nsgfrmmat(g,shift,M,Ls,steps)
%NSGFRMMAT  Sparse nonstationary Gabor frame operator matrix
%   Usage:  S = nsgfrmmat(g,shift,M,Ls)
%           S = nsgfrmmat(g,shift,M)
%           S = nsgfrmmat(g,shift)
%
%   Input parameters:
%         g         : Cell array of window functions
%         shift     : Vector of shifts between the center frequencies
%         M         : Vector of lengths of the window functions
%         Ls        : Signal length (optional)
%         steps     : Maximum number of side-diagonals to compute
%                     (optional)
%   Output parameters:
%         S         : Sparse frame operator matrix
%
%   Creates the frame operator matrix of size LsxLs associated to 
%   the nonstationary Gabor system g, shift, M in sparse matrix 
%   format.
%
%   From the Walnut representation of the nonstationary Gabor frame 
%   operator S we can deduce that the discrete nonstationary Gabor frame 
%   operator is represented by a sparse matrix. 
%   
%   Let N = numel(shift) and K_l = {nin [0,N-1] : l=0 mod M(n)}, 
%   then
%
%       S(k,j) =     sum       M(n) g{n}[k]*conj(g{n}[j]).
%                n in K_{k-j}
%
%   This representation is used together with the size (support) of the 
%   windows g to compute only the relevant entries of S are computed.
%   
%   The optional parameter steps can be used to compute approximations of
%   the frame operator with only a certain number of side-diagonals.
%
%   See also: nsgt, nsigt, nsganamat
%
%   References:
%     N. Holighaus. Structure of nonstationary Gabor frames and their dual
%     systems. Preprint, 2013.
%     
%     P. Balazs, M. Dörfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
%     Theory, implementation and applications of nonstationary Gabor Frames.
%     J. Comput. Appl. Math., 236(6):1481-1496, 2011.
%     
%
%   Url: http://nsg.sourceforge.net/doc/operators/nsgfrmmat.php

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
% Date: 24.04.13

if nargin < 2
    error('Too few input arguments');
end

N = length(g);
posit = cumsum(shift)-shift(1);

Lg = cellfun(@length,g);
g = cellfun(@fftshift,g,'UniformOutput',0);

if nargin < 5
    steps = Inf;
    if nargin < 4
        Ls = posit(end)+shift(1);
        if nargin < 3
            M = Lg;
        end
    end
end

B = floor((Lg-1)./M);
S0 = sparse(Ls,Ls);

for ll = 1:N
    
    win_range = mod(posit(ll)+(-floor(Lg(ll)/2):ceil(Lg(ll)/2)-1),Ls)+1;
    S0(win_range,1) = S0(win_range,1) + M(ll)*abs(g{ll}).^2;
    
    for kk = 1:min(B(ll),steps)
        temp0 = M(ll)*g{ll}(1+kk*M(ll):end).*conj(g{ll}(1:end-kk*M(ll)));
        win_begin = win_range(1:end-kk*M(ll));
        win_end = win_range(1+kk*M(ll):end);
        
        S0(win_begin,1+kk*M(ll)) = S0(win_begin,1+kk*M(ll)) + conj(temp0);
        S0(win_end,Ls+1-kk*M(ll)) = S0(win_end,Ls+1-kk*M(ll)) + temp0;
    end
end

[i,j] = find(S0);
ind0 = sub2ind([Ls,Ls],i,j);
S = sparse(i,mod(i+j-2,Ls)+1,S0(ind0),Ls,Ls);
