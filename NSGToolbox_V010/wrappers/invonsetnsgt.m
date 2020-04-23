function fr = invonsetnsgt(c,g,shift,M,Ls)
%INVONSETNSGT  Onset-based nonstationary Gabor synthesis
%   Usage:  fr = invonsetnsgt(c,g,shift,M,Ls)
%           fr = invonsetnsgt(c,g,shift,M)  
%
%   Input parameters: 
%         c         : Cell array of transform coefficients
%         g         : Cell array of analysis windows
%         shift     : Vector of time shifts
%         M         : Number of time channels
%         Ls        : Original signal length
%   Output parameters: 
%         fr        : Reconstructed signal
%
%   This is a wrapper function for the inverse scaleframe nonstationary 
%   Gabor transform with onset detection based adaptation. It basically
%   just forwards the input to NSDUAL and NSIGT_REAL. 
%
%   For more information see ONSETNSGT.
%
%   See also:  onsetnsgt, nsigt_real, nsdual
%
%   References:
%     P. Balazs, M. Dörfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
%     Theory, implementation and applications of nonstationary Gabor Frames.
%     J. Comput. Appl. Math., 236(6):1481-1496, 2011.
%     
%
%   Url: http://nsg.sourceforge.net/doc/wrappers/invonsetnsgt.php

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
    if nargin < 4
        error('Not enough input arguments');
    end
    Ls = sum(shift);
end

gd = nsdual(g,shift,M);

fr = nsigt_real(c,gd,shift,M,Ls);
