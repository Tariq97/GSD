function [Y,spar,spec,ierr,Rseed] = SpecSyn3(N,samp,corr,acf,Rseed,lmin,lmax)
%
%
% Modified by Tariq  - added positive slip at centre code of Martin Mai
%
%
%  [Y,spar,spec,ierr] = SpecSyn2(N,samp,corr,'acf',Rseed) 
%  
%  generates a 2D-random field Y of size (Nz+1 x Nx+1) with
%  possibly variable spatial sampling in z,x-direction. 
%  This function simulates anisotropic random fields, i.e. 
%  the correlation length in both directions can be different; 
%  rectangular grid dimensions are also possible.
%  corr contains the corr. length ax, az and the Hurstnumber H 
%  or the fractal dimension D; the autocorrelation function to 
%  be used has to be specified by acf;
%  NOTE: D = 3-H, and larger D (smaller H) yield "rougher" fields
%
%  The algorithm is based on the spectral synthesis method by 
%	Pardo-Iguzquiza, E. and Chica-Olma, M. (1993)
%	The Fourier integral method: and efficient spectral method for 
%	simulation of random fields, Mathematical Geology, 25, p177-217.
%  but extends their method to handle rectangular grids and variable
%  sampling in the two directions.
%
%  INPUT:
%  N     - grid dimensions [Nz Nx]
%  samp	 - desired sampling, [dz dx]; if scalar, then dz = dx
%  corr  - corr = [az ax]   for 'gs' and 'ex' 
%	   corr = [az ax H] for 'ak'; note that 0 <= H <= 1
%	   corr = [D kc]    for 'fr'; D is the fractal dimension,
%	    kc: corner wavenumber, spectrum decays linearly for k>kc
%  acf   - autocorrelation function: 
%	  'gs' or 'GS' - Gaussian
%	  'ex' or 'EX' - Exponential
%	  'ak' or 'AK' - anisotropic vonKarman
%	  'fr' or 'FR' - fractal distribution
%  Rseed - seeds for random number generators; if omitted or empty
%	   Rseed = sum(100*clock) is used (returned in structure spar)
%	   [Rpseed Rsseed] for the phase and small random spectral part
% 
%  OUTPUT:
%  Y    - Y = [z x] random field whose size is determined by the sample
%	  spacing in each direction, the number of points and whether
%	  the pow2-option is given or not. 
%  spar	- structure with length vectors, sampling in both direction
%	  and other relevant parameters; it also contains the random
%	  seed number, useful to reproduce realizations
%  spec - structure containing the computed power spectrum as well
%         wavenumber vectors
%  ierr - 0 when successfully executed; 
%         1 when error in Z,X sampling;
%  Rseed (added by Luis)
%
%  Written by Martin Mai (martin@seismo.ifg.ethz.ch) 
%  originally from 07/16/98, based on SRB-toolbox (Ph. Rio)
%  last changes 03/01/2000; Nov. 2002;
%               May 2015, Martin Galis - see *** MG *** below
%
% -------------------------------------------------------------------------

ierr  = 0;    % error variable
check = 'y';  % set to 'y' if you want to create
              % a simple out put plot to check
              % the spectra and the resulting field

% check input variables
if     nargin < 4;  error('Not enough input arguments');
elseif nargin == 4; Rseed = [];
end;

if length(samp) == 1; samp = [samp samp]; end
if length(N) == 1; N = [N N]; end


% error checking on inpur array size and given sampling
if mod(N(2),samp(2)) ~= 0
 ierr = 1;
 disp('** sampling in X does not yield an integer number');
 disp('   of grid points');
 disp('==> BOOM OUT in SpecSyn2<=='); return
end
if mod(N(1),samp(1)) ~= 0
 ierr = 1;
 disp('** sampling in Z does not yield an integer number');
 disp('   of grid points ==> abort!');
 disp('==> BOOM OUT in SpecSyn2<=='); return
end


% get data values on the correlation length/fractal dimension
if  (acf == 'fr') | (acf == 'FR')
  if (length(corr) == 2)
	D = corr(1); kc = corr(2);
  else
	D = corr(1); kc = 0.1; 
	disp('** Corner wavenumber kc not given: set to 0.1 **')
  end
elseif (acf == 'ex') | (acf == 'EX') | (acf == 'gs' ) | (acf == 'GS')
	ax = corr(2);  az = corr(1);
elseif (acf == 'ak') | (acf == 'AK')
	ax = corr(2);  az = corr(1); H = corr(3);
end


% % % set size for spectral synthesis tool that generates
% % % fields of size (2*rmz+1, 2*rmx+1), i.e. the method 
% % % requires an ODD number of points in each direction
%nptsX = N(2)/samp(2);  %% number of grid-points in X
%nptsZ = N(1)/samp(1);  %% number of grid-points in Z
nptsX = round(N(2)/samp(2));  %% number of grid-points in X
nptsZ = round(N(1)/samp(1));  %% number of grid-points in Z
%if     mod(nptsX,2) == 0, rmx = nptsX/2;
%elseif mod(nptsX,2) == 1, rmx = (nptsX-1)/2;
%end
%if     mod(nptsZ,2) == 0, rmz = nptsZ/2;
%elseif mod(nptsZ,2) == 1, rmz = (nptsZ-1)/2;
%end

if     mod(nptsX,2) == 0, rmx = nptsX/2;
elseif mod(nptsX,2) ~= 0, rmx = (nptsX-1)/2;
end
if     mod(nptsZ,2) == 0, rmz = nptsZ/2;
elseif mod(nptsZ,2) ~= 0, rmz = (nptsZ-1)/2;
end

% % % compose power spectrum for two of the four quadrants
% % % wavenumber vector in [-pi,pi]
kx = (rmx./((2*rmx+1)*samp(2)))*linspace(-2*pi,2*pi,2*rmx+1);
kz = (rmz./((2*rmz+1)*samp(1)))*linspace(-2*pi,2*pi,2*rmz+1);
kr = zeros(rmz+1,2*rmx+1);

for j = 1:(2*rmx+1)
   for i = 1:(rmz+1)
	if ((acf == 'fr') | (acf == 'FR'))
	   kr(i,j) = sqrt((kz(i)^2) + (kx(j)^2));
	else
	   kr(i,j) = sqrt((az^2)*(kz(i)^2) + (ax^2)*(kx(j)^2));
	end
   end
end


% % % calculate power spectral density, depending on selected ACF
if ((acf == 'gs') | (acf == 'GS')) 
	PS = 0.5*ax*az*exp(-0.25*kr.^2);
elseif ((acf == 'ex') | (acf == 'EX'))
 	PS = (ax * az)./(1 + (kr.^2)).^1.5;
elseif ((acf == 'ak') | (acf == 'AK'))
	k = kr(:); k = k(k>0);
	coef = 4*pi*H*ax*az./besselk(H,min(k));
%	coef = ax*az; 
  	PS = coef./(1 + (kr.^2)).^(H+1);
elseif ((acf == 'fr') | (acf == 'FR'))
	decay = 0.5*(8-2*D);
	% % to ensure proper scaling of the power spectrum 
	% % in k-space we cannot allow kr == 0
	if min(kr(:)) == 0
	  [p,q] = find(kr == 0);
	  k(p,q) = mean(mean(kr(p-1:p,q-1:q)));
    end
	% % set values below k< kc to constant pi*kc
	%kr(kr <= kc) = pi*kc;		
	% % set values below k< kc to constant pi*kc -- this caused a peak in
	% PS, value kc is replaced by pi*kc...  
	kr(kr <= kc) = kc;		

    PS = 1./((kr.^2).^decay);         	
end

% remove wavelength shorter than lmin, i.e., k larger than 2*pi/lmin:
if lmin > 0
    kmax = 2.0*pi/lmin;
    PS(:,abs(kx)>kmax) = 0;
    PS(kz(1:rmz+1)<-kmax,:) = 0;
end

% % % remove wavelength longer than lmax, i.e., k smaller than 2*pi/lmax:
% % if lmax > 0
% %     kmin = 2.0*pi/lmax;
% %     PS(:,abs(kx)<kmin) = 0;
% %     figure; loglog(kx,PS(end,:));
% %     PS(kz(1:rmz+1)>-kmin,:) = 0;
% %     figure; loglog(kz(1:rmz+1),PS(:,601));
% %     
% % %     PS(:,:) = 0;
% %     
% %     PS(rmz+1,600) = 1e8;
% %     PS(rmz+1,602) = 1e8;
% % 
% %     PS(300,601) = 1e8;
% % 
% % end


% % % the IFFT needs the spectrum normalized to max.amp unity
PS = PS./max(PS(:));
AM = sqrt(PS);

% % % compose the random phase
% % % initialize random number generator
if isempty(Rseed) == 1;
 %Rseed = zeros(1,2);
 Rseed(1) = sum(100*clock);
 %Rseed(2) = sum(109*clock);
end
rand('seed',Rseed(1));
%randn('seed',Rseed(2));

% % % random phase in [0,pi]
PH = 2*pi*rand(size(kr));	 

% % % assemble the random field in FT-domain	
% add small random high-wavenumber components
%x = (1 + 0.5*randn(size(kr)));
x = 1;
RAD = AM.*x;

% % % set DC-component to different value, if desired
% % % NOTE that this only changes the overall 'level' of
% % % the field, not its appearance, but has significant
% % % impact on the Fourier transform which reflects the
% % % DC-value at the smallest wavenumber ("Nugget Effect")
Neff = 0;                               % % "Nugget" alue
RAD(rmz+1,rmx+1) = Neff;	          	% % "Nugget" effect, zero-mean field
AM(rmz+1,rmx+1)  = RAD(rmz+1,rmx+1);
Y = RAD.*cos(PH) + sqrt(-1)*RAD.*sin(PH);

% % % the following takes care of the conjugate symmetry condition
% % % in order to ensure the final random field is purely real
U = zeros(2*rmz+1,2*rmx+1);	        % % will be conj. sym. field
Y = [Y ; conj(fliplr(flipud(Y(1:rmz,:))))];

for i = 1:rmx
    Y(rmz+1,-i+2*rmx+2) = conj(Y(rmz+1,i));
end
for i = 1:rmz+1
  for j = 1:rmx+1
     U(i,j) = Y(i+rmz,j+rmx);
   end
end
for i = rmz+2:2*rmz+1
  for j = rmx+2:2*rmx+1
     U(i,j) = Y(i-rmz-1,j-rmx-1);
   end
end
for i = 1:rmz+1
  for j = rmx+2:2*rmx+1
     U(i,j) = Y(i+rmz,j-rmx-1);
   end
end
for i = rmz+2:2*rmz+1
  for j = 1:rmx+1
     U(i,j) = Y(i-1-rmz,j+rmx);
   end
end

% % % take 2D-inverse FFT to obtain spatial field; imaginary parts
% % % of order 1e-13 due to machine precision are removed; 
% % % also, remove mean and scale to unit variance
Y = real(ifft2(U));
Y = Y/std2(Y);				% standard deviation of unity


[lz,lx] = size(Y); 
px = round(lx/sqrt(2));		% dimensions for 'interior' area
pz = round(lz/sqrt(2));	
qx = floor(0.5*(lx-px));	% indices for 'interior' area
qz = floor(0.5*(lz-pz));		
GI = Y(qz+1:end-qz,qx+1:end-qx);
if mean(GI(:)) < mean(Y(:))
  Y = -1*Y; 
end

spar.dim   = N;
spar.samp  = samp; 
spar.size  = size(Y);
spar.corr  = corr;
spar.acf   = acf;
spar.Rseed = Rseed;
spar.lx = [0:samp(2):samp(2)*(size(Y,2)-1)]; 
spar.lz = [0:samp(1):samp(1)*(size(Y,1)-1)];

px = find(kx >= -1e-8);
pz = find(kz <=  1e-8);
spec.PD  = PS(:,px);
spec.kpx = kx(px);
spec.kpz = kz(pz);
spec.PDx = spec.PD(end-1,:);
spec.PDz = spec.PD(:,1);


