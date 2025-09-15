function [S] = TaperSlip(S,N,window)
%
% function [S] = TaperSlip(S,N,'window')
% tapers the slip amplitudes at the fault boundaries 
% to avoid large slip at the edges. 
%	
% INPUT:
% S	- original grid
% N 	- vector of numbers of rows/cols to be tapered [left/right top bottom]
% window- 'hn' for Hanning window 
%	  'kw' for Kaiser window
%	  'tr' for triangular window
% 	  window can also be a vector [c1 c2 c3] with c(i) < 1; 
%	  in this case, the corresponding rows/cols will simply be 
%	  multiplied by that c(i)
%	  NOTE: this scheme only works if only one col/row is tapered
%
% Default values: N = [1 0 1], window = 'hn'
%
% OUTPUT: 
% S 	- tapered grid (same size as input grid)
%
% See also MATLAB function HANNING, and KAISER
	
% Written by Martin Mai (mmai@pangea.Stanford.EDU) 
% 05/08/98
% last change 08/26/99
% ------------------------------------------------

if nargin == 1;
	N = [1 0 1];
	window = 'hn';
elseif nargin == 2;	
	window = 'hn';
end

%%% create taper window, i.e. Kaiser window (NOTE: increasing beta in 
%%% Kaiser-Window widens main lobe and attenuates amplitude in the side lobes)
%%% or hanning (cosine) from 0 to 1 in N steps on EACH side
if isstr(window) == 1
  if window == 'hn'
	taperS = hanning(2*N(1)+1);	% for left/right columns
	taperT = hanning(2*N(2)+1);	% for top row
	taperB = hanning(2*N(3)+1);	% for bottom rows
  elseif window == 'kw'
	beta = 6;
	taperS = kaiser(2*N(1)+1,beta);
	taperT = kaiser(2*N(2)+1,beta);
	taperB = kaiser(2*N(3)+1,beta);
  elseif window == 'tr'
	taperS = triang(2*N(1)+1);	
	taperT = triang(2*N(2)+1);	
	taperB = triang(2*N(3)+1);	
  end
  winS = taperS(N(1)+2:2*N(1)+1)'; 
  winT = taperT(1:N(2))'; 
  winB = taperB(end-(N(3)-1):end)'; 

elseif isstr(window) == 0
    [i,j] = find(N == 0);  % to make sure that rows/cols with N = 0 are not
  window(j) = 1;	   % tapered in case s contains entries other than 1	
  winS = window(1);
  winT = window(2);
  winB = window(3);

end

%%% check in case no taper is applied at one of the boundaries
if isempty(winS) == 1; winS = 1; end	
if isempty(winT) == 1; winT = 1; end
if isempty(winB) == 1; winB = 1; end

 
%%% set-up bell-shaped tapering function 
bell = ones(size(S)); 
[j,k] = size(S);
ls = length(winS); 
lt = length(winT);
lb = length(winB);
for zs = 1:j;  
    bell(zs,1:ls) = bell(zs,1:ls).*fliplr(winS); 
    bell(zs,k-ls+1:k) = bell(zs,k-ls+1:k).*winS;
end
for zt = 1:k
    bell(1:lt,zt) = bell(1:lt,zt).*winT(:);
    bell(j-lb+1:j,zt) = bell(j-lb+1:j,zt).*winB(:);
end


%%% finally, multiply input slip with bell-function
S = S.*bell;
