function fr = nsigt_real(c,g,shift,M,Ls)
% NSIGT_REAL  Nonstationary Gabor synthesis for real signals
%   Usage: fr = nsigt_real(c,g,shift,M,Ls)
%
%   Input parameters: 
%         c         : Cell array of nonstationary Gabor coefficients
%         g         : Cell array of synthesis windows
%         shift     : Vector of time shifts
%         M         : Number of frequency channels (vector/scalar)
%         Ls        : Length of the analyzed signal
%   Output parameters:
%         fr        : Synthesized real-valued signal (Channels are stored 
%                     in the columns)
%
%   Given the cell array c of nonstationary Gabor coefficients, a set of 
%   windows g and time shifts shift, this function computes the 
%   corresponding real-valued nonstationary Gabor synthesis. Let 
%   N=numel(g) and P(n)=sum_{l=1}^{n} shift(l), then the complex 
%   valued synthesis formula reads:
%
%               N-1 
%       fr(l) = sum sum c{n}(m)g{n}[l-P(n)]*exp(2*pi*i*(l-P(n))*m/M(n)),
%               n=0  m
%   
%   for l=0,cdots,Ls-1. In practice, the synthesis formula is realized 
%   by ifft and overlap-add. In the real valued case, ifftreal provides
%   the missing frequency content normally given by the coefficients
%   c{n}(m) for m >= floor(M(n)/2).
% 
%   If a nonstationary Gabor frame was used to produce the coefficients 
%   and g is a corresponding dual frame, this function should perfectly 
%   reconstruct the originally analyzed signal to numerical precision.
%
%   Note that NSIGT_REAL requires the input parameter M to guarantee
%   that the vectors used in the overlap-add process are of the correct
%   length.
%
%   Multichannel output will save each channel in a column of fr.
%
%   See also:  nsgt_real, nsdual, nstight
% 
%   References:
%     P. Balazs, M. Dörfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
%     Theory, implementation and applications of nonstationary Gabor Frames.
%     J. Comput. Appl. Math., 236(6):1481-1496, 2011.
%     
%
%   Url: http://nsg.sourceforge.net/doc/core_routines/nsigt_real.php

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

% Author: Nicki Holighaus, Gino Velasco
% Date: 23.04.13

% Preparation

if nargin < 4
    error('Not enough input arguments');
end

if iscell(c) == 0 % If matrix format coefficients were used, convert to
    % cell
    [~,N,CH] = size(c);
    c = reshape(c,N*M,CH);
    c = mat2cell(c,M*ones(N,1),CH);
    M = M*ones(N,1);
else
    N = length(c);
    CH = size(c{1},2);
end

posit = cumsum(shift);        % Calculate positions from shift vector
NN = posit(end);              % Reconstruction length before truncation
posit = posit-shift(1);   	% Adjust positions

fr = zeros(NN,CH); % Initialize output

if nargin < 5
    Ls = NN; % If original signal length is not given do not truncate
end

% The overlap-add procedure including multiplication with the synthesis
% windows

for ii = 1:N
    Lg = length(g{ii});
    
    win_range = mod(posit(ii)+(-floor(Lg/2):ceil(Lg/2)-1),NN)+1;
    
    temp = ifftreal(c{ii},M(ii),1)*M(ii);
    temp = temp(mod([end-floor(Lg/2)+1:end,1:ceil(Lg/2)]-1,...
        length(temp))+1,:);
    
    fr(win_range,:) = fr(win_range,:) + ...
        bsxfun(@times,temp,g{ii}([Lg-floor(Lg/2)+1:Lg,1:ceil(Lg/2)]));
end

fr = fr(1:Ls,:); % Truncate the signal to original length (if given)
