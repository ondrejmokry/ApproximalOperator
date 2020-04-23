function [c,Ls] = nsgt(f,g,shift,M)
%NSGT  Nonstationary Gabor transform
%   Usage: [c,Ls] = nsgt(f,g,shift,M)
%          [c,Ls] = nsgt(f,g,shift)
%          c = nsgt(...)
%
%   Input parameters:
%         f         : The signal to be analyzed (For multichannel
%                     signals, input should be a matrix which each
%                     column storing a channel of the signal).
%         g         : Cell array of frequency side analysis windows
%         shift     : Vector of time shifts
%         M         : Number of frequency channels (optional).
%                     If M is constant, the output is converted to a
%                     matrix
%   Output parameters:
%         c         : Transform coefficients (matrix or cell array)
%         Ls        : Original signal length (in samples)
%
%   Given the cell array g of windows, the time shift vector shift, and
%   channel numbers M, NSGT computes the corresponding nonstationary 
%   Gabor transform of f. Let P(n)=sum_{l=1}^{n} shift(l), then the 
%   output c = NSGT(f,g,shift,M) is a cell array with 
%
%              Ls-1                               
%      c{n}(m)= sum f(l)*conj(g\{n\}(l-P(n)))*exp(-2*pi*i*(l-P(n))*m/M(n))
%               l=0                                
%
%   where m runs from 0 to M(n)-1.
%
%   If multichannel input is used, the same nonstationary Gabor system is 
%   applied to each channel and each entry of c will be a 2D array with 
%   c{n}(:,CH) being the entries corresponding to time position n, 
%   signal channel CH.
%
%   If M is scalar or uniform, then c is converted into a regular
%   array. 
%
%   The choice of phase-locked coefficients (by inserting (l-P(n)) in the 
%   complex exponential) prevents border artifacts for combinations of 
%   window functions g{n} that do not have full support and M(n) that 
%   do not divide Ls.
%
%   See also:  nsigt, nsdual, nstight
%
%   References:
%     P. Balazs, M. Dörfler, F. Jaillet, N. Holighaus, and G. A. Velasco.
%     Theory, implementation and applications of nonstationary Gabor Frames.
%     J. Comput. Appl. Math., 236(6):1481-1496, 2011.
%     
%
%   Url: http://nsg.sourceforge.net/doc/core_routines/nsgt.php

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

% Check input arguments

if nargin <= 2
    error('Not enough input arguments');
end

[Ls,CH]=size(f);

if Ls == 1
    f = f.';
    Ls = CH;
    CH = 1;
end

if CH > Ls
    disp(['The number of signal channels (',num2str(CH),') ',...
        'is larger than']);
    disp(['the number of samples per channel (',num2str(Ls),').']);
    reply = input('Is this correct? ([Y]es,[N]o)', 's');
    switch reply
        case {'N','n','No','no',''}
            reply2 = input('Transpose signal matrix? ([Y]es,[N]o)', 's');
            switch reply2
                case {'N','n','No','no',''}
                    error('Invalid signal input, terminating program');
                case {'Y','y','Yes','yes'}
                    disp('Transposing signal matrix and continuing ',...
                        'program execution');
                    f = f.';
                    X = CH; CH = Ls; Ls = CH; clear X;
                otherwise
                    error('Invalid reply, terminating program');
            end
        case {'Y','y','Yes','yes'}
            disp('Continuing program execution');
        otherwise
            error('Invalid reply, terminating program');
    end
end

N=length(shift);    % The number of time slices

if nargin == 3
    M = zeros(N,1);
    for kk = 1:N
        M(kk) = length(g{kk});
    end
end

if max(size(M)) ==  1
    M = M(1)*ones(N,1);
end

posit = cumsum(shift)-shift(1); % Calculate positions from shift vector

% A small amount of zero-padding might be needed (e.g. for scale frames)

fill = sum(shift)-Ls;
f = [f;zeros(fill,CH)];

c=cell(N,1); % Initialisation of the result

% The actual transform

for ii = 1:N
    Lg = length(g{ii});
    
    idx = [ceil(Lg/2)+1:Lg,1:ceil(Lg/2)];
    win_range = mod(posit(ii)+(-floor(Lg/2):ceil(Lg/2)-1),Ls+fill)+1;
    
    if M(ii) < Lg % if the number of frequency channels is too small,
        % aliasing is introduced
        col = ceil(Lg/M(ii));
        temp = zeros(col*M(ii),CH);
        
        temp([end-floor(Lg/2)+1:end,1:ceil(Lg/2)],:) = ...
            bsxfun(@times,f(win_range,:),g{ii}(idx));
        temp = reshape(temp,M(ii),col,CH);
        
        c{ii} = squeeze(fft(sum(temp,2)));
        % Using c = cellfun(@(x) squeeze(fft(x)),c,'UniformOutput',0);
        % outside the loop instead does not provide speedup; instead it is
        % slower in most cases.
    else
        temp = zeros(M(ii),CH);
        temp([end-floor(Lg/2)+1:end,1:ceil(Lg/2)],:) = ...
            bsxfun(@times,f(win_range,:),g{ii}(idx));
        
        c{ii} = fft(temp);
    end
end

if max(M) == min(M)
    c = cell2mat(c);
    c = reshape(c,M(1),N,CH);
end
