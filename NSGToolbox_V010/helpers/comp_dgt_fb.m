function [coef]=comp_dgt_fb(f,g,a,M)
%COMP_DGT_FB  Filter bank DGT
%   Usage:  c=comp_dgt_fb(f,g,a,M);
%  
%   This is a computational routine. Do not call it directly.
%
%   This code was borrowed from LTFAT with kind permission of author 
%   Peter L. Søndergaard. 
%
%   Please visit the LTFAT toolbox homepage:
%   http://ltfat.sourceforge.net/
%
%   Url: http://nsg.sourceforge.net/doc/helpers/comp_dgt_fb.php

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

%   AUTHOR : Peter L. Søndergaard.
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version. To view the license, see 
%   <http://www.gnu.org/licenses/>.


% Calculate the parameters that was not specified.
L=size(f,1);
N=L/a;
gl=length(g);
W=size(f,2);      % Number of columns to apply the transform to.
glh=floor(gl/2);  % gl-half


% Conjugate the window here.
g=conj(fftshift(g));

coef=zeros(M,N,W);

% ----- Handle the first boundary using periodic boundary conditions. ---
for n=0:ceil(glh/a)-1

    % Periodic boundary condition.
    fpart=[f(L-(glh-n*a)+1:L,:);...
           f(1:gl-(glh-n*a),:)];
    
    fg=bsxfun(@times,fpart,g);
    
    % Do the sum (decimation in frequency, Poisson summation)
    coef(:,n+1,:)=sum(reshape(fg,M,gl/M,W),2);
      
end;

% ----- Handle the middle case. ---------------------
for n=ceil(glh/a):floor((L-ceil(gl/2))/a)
  
  fg=bsxfun(@times,f(n*a-glh+1:n*a-glh+gl,:),g);
  
  % Do the sum (decimation in frequency, Poisson summation)
  coef(:,n+1,:)=sum(reshape(fg,M,gl/M,W),2);
end;

% ----- Handle the last boundary using periodic boundary conditions. ---
for n=floor((L-ceil(gl/2))/a)+1:N-1

    % Periodic boundary condition.
    fpart=[f((n*a-glh)+1:L,:);... %   L-n*a+glh elements
           f(1:n*a-glh+gl-L,:)];  %  gl-L+n*a-glh elements      
    
    fg=bsxfun(@times,fpart,g);
    
    % Do the sum (decimation in frequency, Poisson summation)
    coef(:,n+1,:)=sum(reshape(fg,M,gl/M,W),2);      
end;

% --- Shift back again to make it a frequency-invariant system. ---
for n=0:N-1
  coef(:,n+1,:)=circshift(coef(:,n+1,:),n*a-glh);
end;

coef=fft(coef);
