function peaks = peakpick(SF,thre,range,multi)
%PEAKPICK  Peakpicking routine
%   Usage:  peaks = peakpick(SF,thre,range,multi)
%
%   Input parameters: 
%         SF        : Onset detection function
%         thre      : Threshold value
%         range     : Relevance area for local maximum
%         multi     : Asymmetric extension factor for relevance area
%   Output parameters:
%         peaks     : Significant maxima of SF*
%
%   This is a helper function for ONSETDET and not meant to be used 
%   individually.
%
%   For an onset detection function SF, the routine picks only those 
%   local maxima that are larger than the local mean over an area of the 
%   form
%
%           <---multi*range---X---range--->
%
%   by more than the threshold given by thre.
%
%   See also:  onsetdet
%
%   References:
%     S. Dixon. Onset detection revisited. In Proceedings of the 9th
%     International Conference on Digital Audio Effects, volume 120, pages
%     133-137, 2006.
%     
%
%   Url: http://nsg.sourceforge.net/doc/helpers/peakpick.php

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

if nargin < 4
    error('Not enough input arguments');
end

% Compute local maxima

maxima = locmax(SF,range);

m = length(maxima);
n = length(SF);

% Since the signal is considered periodic, some
% periodic border values are added

SF = [SF(n-(multi+1)*range+1:1:n),SF,SF(1:range)];

kk=0;
peaks = [];

% The loop selects the significant local maxima
% from the output of locmax

for ii = 1:m
    pos = maxima(ii);
    th_loc_mean = sum(SF(pos:pos+(multi+2)*range))/((multi+2)*range+1)+thre;
    if ( SF(pos+(multi+1)*range) > th_loc_mean )
        kk = kk + 1;
        peaks(kk) = maxima(ii);
    end
end   

end

% Computes 'some kind of' local maxima.
% A point of SF is considered a local maximum
% in the sense of this routine, if it is larger
% than all surrounding points in a neigborhood
% with a radius of range-1 of this point.

function maxima = locmax(SF,range)

if ( size(SF,1) ~= 1 )
    SF=SF.';
end

n = length(SF);

SF = [SF(n-range+1:1:n),SF,SF(1:range)];

kk = 1;
maxima = [];

for ii = 1:n
    xx = SF(ii:ii+2*range);
    yy = ones(1,2*range+1)*SF(ii+range);
    if ( sum( yy >= xx ) == 2*range+1 )
        maxima(kk) = ii; 
        kk = kk+1;
    end   
end

end
