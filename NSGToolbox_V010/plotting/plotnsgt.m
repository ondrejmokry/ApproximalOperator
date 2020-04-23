function plotnsgt(c,shift,sr,varargin)
%PLOTNSGT  Plot nonstationary Gabor coefficients
%   Usage:  plotnsgt(c,shift,sr,varargin)
%           plotnsgt(c,shift,sr)
%           plotnsgt(c,shift)
%
%   Input parameters:
%         c        : Array of coefficients.
%         shift    : Vector of time shifts
%         sr       : Sampling rate in Hz (default 1 Hz)
%         varargin : Optional input parameters (see table below)
%
%   Given a coefficient array c and the time shift vector shift, this
%   function plots the dB-scaled nonstationary Gabor spectrogram 
%   corresponding to c. To capture the correct position of the 
%   coefficients in the time frequency plane, the columns of the 
%   spectrogram (coefficients corresponding to the same time position) are 
%   stretched accordingly.
%
%   If additionally, the sampling rate sr is provided, time and frequency 
%   axes will be labeled properly. 
%
%   If the coefficients were obtained using NSGT_REAL, the realsig*
%   switch should be used, otherwise only half the desired frequency range
%   will be displayed. The shown frequency range can be controlled with the
%   cutout parameter (default: 2) and the dynamic range of the
%   spectrogram can be adjusted with dynrange.
%
%   Optional input arguments arguments can be supplied like this:
%
%       plotnsgt(c,shift,sr,'dynrange',dynrange)
%
%   The arguments must be character strings followed by an
%   argument:
%
%     'dynrange',dynrange  Colorscale dynamic range in dB (default 60 dB)
%
%     'cutout',cutout      Desired part of the spectrogram, e.g.
%                          choice of '2' shows frequencies up to Nyquist
%                          ('X' shows the 'number_of_bins/X' lowest 
%                          frequency bins)
%
%     'realsig',realsig    Input coefficients are taken from a
%                          representation for real-valued signals
%
%   See also:  nsgt, plotnsgt
%
%   Url: http://nsg.sourceforge.net/doc/plotting/plotnsgt.php

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

% Author:  Nicki Holighaus
% Original code by Florent Jaillet
% Date: 26.04.13

if nargin < 3
    % Default value for sampling frequency.
    sr = 1;
    if nargin < 2
        error('Not enough input arguments');
    end
end

dynrange = 60;  % Default value for colorscale dynamic.
cutout = 2;     % Default value for frequency cutout.
realsig = 0;

if nargin >= 4
    Lvar = length(varargin);
    if mod(Lvar,2)
        error('Invalid input argument');
    end
    for kk = 1:2:Lvar
        if ~ischar(varargin{kk})
            error('Invalid input argument');
        end
        switch varargin{kk}
            case {'dynrange'}
                dynrange = varargin{kk+1};
            case {'cutout'}
                cutout = varargin{kk+1};
            case {'realsig'}
                realsig = varargin{kk+1};
            otherwise
                error(['Invalid input argument: ', varargin{kk}]);
        end
    end
end

maxfreq = sr/cutout;

if realsig
    cutout = cutout/2;
end

cla %Clear previous axes

posit=cumsum(shift)-shift(1);

% Compute time limit for the representation of each window.
tlim=diff(posit)/2;
tlim=[posit(1)-tlim(1);posit(1:end-1)+tlim;posit(end)+tlim(end)];
tlim=tlim/sr;

% Compute maximum of the representation for colorscale dynamic handling.
if iscell(c) == 1
    if size(c{1},2) > 1
        error(['Multichannel spectrograms are not supported. Please ',...
            'use ''cellfun(@(x) x(:,k),c,''UniformOutput'',0)'' to ',...
            'select the k-th channel.']);
    end
    temp=cell2mat(c);
else
    if size(c,3) > 1
        error(['Multichannel spectrograms are not supported. Please ',...
            'use ''c(:,:,k)'' to select the k-th channel.']);
    end
    temp=c;
end
ma=20*log10(max(abs(temp(:))));

% Plot the representation: as the sampling grid in the time frequency plane
% is irregular, the representation by done by plotting many images next to
% each other, with one image for each window
hold('on');
if iscell(c) == 1
    for ii=1:length(shift)
        temp = 20*log10(abs(c{ii})+eps);
        ind = ceil(0.5:0.5:length(temp)/cutout);
        temp = temp(ind);% +eps is here to avoid log of 0
        % Octave cannot plot images that are only one point wide, so we use
        % images that are to points wide
        imagesc(adapt(tlim(ii:ii+1)),[0,1-1/length(c{ii})]*maxfreq,...
            [temp,temp],[ma-dynrange,ma]);
    end
else
    for ii=1:length(shift)
        temp = 20*log10(abs(c(ii,:)+eps)).';
        ind = ceil(0.5:0.5:length(temp)/cutout);
        temp = temp(ind);% +eps is here to avoid log of 0
        % Octave cannot plot images that are only one point wide, so we use
        % images that are to points wide
        imagesc(adapt(tlim(ii:ii+1)),[0,1-1/size(c,2)]*maxfreq,...
            [temp,temp],[ma-dynrange,ma]);
    end
end
hold('off');
axis('tight');

end

function [res]=adapt(lim)
% we have to adapt the time values to fit the way the image function handle
% the x position
res=[(3*lim(1)+lim(2))/4,(lim(1)+3*lim(2))/4];
end

