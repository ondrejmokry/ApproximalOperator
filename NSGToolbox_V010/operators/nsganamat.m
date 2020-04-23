function G = nsganamat(g,shift,M,Ls,phaselock)
%NSGANAMAT  Nonstationary Gabor analysis operator matrix
%   Usage: G = nsganamat(g,shift,M,Ls,phaselock)
%          G = nsganamat(g,shift,M,Ls)
%          G = nsganamat(g,shift,M)
%          G = nsganamat(g,shift)
%          G = nsganamat(g,shift,M,phaselock)
%          G = nsganamat(g,shift,phaselock)
%   
%   Input parameters: 
%         g         : Cell array of analysis windows
%         shift     : Vector of time shifts
%         M         : Number of frequency channels (optional)
%         Ls		: Transform length
%         phaselock : This optional 0/1 switch specifies the phaselock 
%                     convention: 0 (non-phaselocked), 
%                     1 (phaselocked, default) 
%   Output parameters:
%         G		    : Frame analysis operator corresponding to the
%                     input arguments
%
%   Computes the frame analysis matrix corresponding the nonstationary
%   Gabor system specified by g, shift and M. The rows of the frame
%   analysis matrix contain the complex conjugate of the frame elements in
%   modulation first order, i.e. let P(n) = sum_{l=1}^{n} shift(l)
%
%             n-1
%       K_n = sum M(l)
%             l=0
%
%       math:: K_n = \sum_{l=0}^{n-1} M(l)
%
%   and let m run from 0 to M(n)-1. Then 
%
%       G(K_n + m,.) = conj(g{n})(. - P(n-1))*exp(-2*pi*i*(. - P(n-1))*m/M(n))
%
%   The conjugate transpose of G equals the synthesis operator
%   corresponding to g, shift and M. Consequently, c=Gf and 
%   fr=conj(G)^T c.
%
%   The formulas above use the phaselocked definition of a nonstationary
%   Gabor frame also realized by NSGT, NSIGT. Alternatively, 
%   non-phaselocked frame elements can be used. Note however, that this
%   might result in border discontinuities if L/M(n) is not integer.
%  
%   Note: While this routine can be used to gain some insight into 
%   the structure of frame-related operators, it is not suited for use with
%   transform lengths over a few thousand samples.
%   
%   See also:  nsgt, nsigt, nsgfrmmat
%
%   References:
%     O. Christensen. Frames and Bases. An Introductory Course. Applied and
%     Numerical Harmonic Analysis. Basel Birkhäuser, 2008.
%     
%     P. Balazs, M. Dörfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
%     Theory, implementation and applications of nonstationary Gabor Frames.
%     J. Comput. Appl. Math., 236(6):1481-1496, 2011.
%     
%
%   Url: http://nsg.sourceforge.net/doc/operators/nsganamat.php

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

% Check input

if nargin == 4
    if Ls == 1;
        phaselock = 1;
        Ls = sum(shift);
    elseif Ls == 0;
        phaselock = 0;
        Ls = sum(shift);
    else
        phaselock = 1;
    end
end
if nargin == 3
    Ls = sum(shift);
    if sum(M) == 1
        phaselock = 1;
        for kk = 1:length(shift)
            M(kk) = length(g{kk});
        end
    elseif sum(M) == 0
        phaselock = 0;
        for kk = 1:length(shift)
            M(kk) = length(g{kk});
        end
    else
        phaselock = 0;
    end
end
if nargin < 3
    if nargin < 2
        error('Not enough input arguments');
    end
    phaselock = 1;
    Ls = sum(shift);
    for kk = 1:length(shift)
        M(kk) = length(g{kk});
    end
end

if max(size(M)) == 1
    M = M(1)*ones(length(shift),1);
end

G = zeros(sum(M),Ls);

% Setup the necessary parameters
N = length(shift);

posit = cumsum(shift);
NN = posit(end);
posit = posit-shift(1);


% Construct the analysis operator matrix explicitly
rows = [1;1+cumsum(M)];
for ii = 1:N
    Lg = length(g{ii});
    
    win_range = mod(posit(ii)+(-floor(Lg/2):ceil(Lg/2)-1),NN)+1;
    G(rows(ii),win_range) = g{ii}([Lg-floor(Lg/2)+1:Lg,1:ceil(Lg/2)]);
    
    if phaselock == 1
        modulates = exp(-2*pi*i*(1:M(ii)-1)'*...
            (-floor(Lg/2):ceil(Lg/2)-1)/M(ii));
        G(rows(ii)+(1:M(ii)-1),win_range) = ...
            bsxfun(@times,modulates,G(rows(ii),win_range));
    elseif phaselock == 0
        modulates = exp(-2*pi*i*(1:M(ii)-1)'*(0:Ls-1)/M(ii));
        G(rows(ii)+(1:M(ii)-1),:) = ...
            bsxfun(@times,modulates,G(rows(ii),:));
    else
        error('Invalid phaselock parameter');
    end
end
