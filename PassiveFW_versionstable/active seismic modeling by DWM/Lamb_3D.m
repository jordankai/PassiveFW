function [T, Urx, Utx, Uzz, Urz] = Lamb_3D(Cs,Cp,r,tmax);
% Point loads suddenly applied onto the surface of a lower elastic halfspace.
% Time variation of load is a unit step function (Heaviside), which is
% applied at the origin on the surface (x=0, z=0).
% Both horizontal and vertical loads are considered.
%
% **********************************************************
%       Version 3.3, September 6, 2012
%       Copyleft by Eduardo Kausel
%       MIT, Room 1-271, Cambridge, MA 02139
%       kausel@mit.edu
% **********************************************************
%
%  Input arguments:
%        pois   = Poisson's ratio
%        NoPlot = plots omitted if true (1). Default is false (0)
%  Output arguments (cylindrical coordinates)
%        T    = Dimensionless time vector, tau = t*Cs/r
%        Urx  = Radial displacement caused by a horizontal load
%               Varies as  cos(theta) with azimuth
%        Utx  = Tangential displacement caused by a horizontal load
%               Varies as  -sin(theta) with azimuth
%        Uzz  = Vertical displacement caused by a vertical load
%               (no variation with azimuth)
%        Urz  = Radial displacement caused by a vertical load
%               (no variation with azimuth)
%        Uzx  = Vertical displacement caused by a horizontal load
%               Not returned, because  Uzx = -Urz (reciprocity),
%               except that Uzx varies as  cos(theta)
%  Note: All displacements are normalized by the shear modulus (mu) and
%        by the range (r), i.e. Urx = urx*r*mu and so forth.
%        Also, the time returned is dimensionless, i.e. T = t*Cs/r, where
%        Cs is the shear wave velocity. Both of these normalizations
%        are equivalent to assuming that mu=1, Cs=1, r=1.
%        Observe that the first point of all arrays returned is for time
%        t=0, and thus not equally spaced with the other points, which
%        begin with the arrival of the P wave. 
% 
% Sign convention:
% Vertical load and vertical displacements at the surfce point up
% Response given in terms of dimensionless time tau=t*Cs/r
%
% References:
%    a) Vertical loads:
%			Eringen and Suhubi, Elastodynamics, Vol. II, 748-750, pois<0.26
% 			Mooney, 1974, BSSA, V. 64, No.2, pp. 473-491, pois > 0.2631,
%           but vertical component only
%    b) Horizontal loads: (pois = 0.25 only)
%        Chao, C.C., Dynamical response of an elastic halfspace to
%        tangential surface loading, Journal of Applied Mechanics,
%        Vol 27, September 1960, pp 559-567
%    c) Generalization of all of the above to arbitrary Poisson's ratio:
%        Kausel, E., Lamb's problem at its simplest
%        Proceeding Royal Society of London, Series A, 2012 (in print)

% Default data & input parameters
pois=getbos(Cs,Cp);
rho = 2000;            % mass density
mu = Cs.^2.*rho;             % shear modulus
             % range (= radial distance)

		% Maximum time for plotting (t=1 => arrival of S waves)
tmax=tmax*Cs/r;
Nt = 10000;          % Number of time steps between tp and ts
NoPlot=1; 


% Roots of Rayleigh function
a2 = (1-2*pois)./(2 - 2*pois);  % (Cs/Cp)^2
b2 = 1-a2;
p = [-16*(1-a2), 8*(3-2*a2), -8, 1];    % Characteristic polynomial
x = sort(roots(p)); % find and sort roots
x1 = x(1);  % First false root
x2 = x(2);  % Second false root
x3 = x(3);  % True Rayleigh root = (Cs/Cr)^2

% Dimensionless arrival times and time arrays
tp = sqrt(a2);		% Arrival time of P waves
ts = 1;				%    "      "   " S   "
tr = sqrt(x3);		%    "      "   " R   "
dt = (ts-tp)/Nt;    % Time step
T0 = [0 tp];        % Time before arrival of P wave
U0 = [0 0];         % Quiescent phase (during T0)
T1 = [tp+dt:dt:ts];     % Time from P to S wave arrival
T2 = [ts+dt:dt:tr-dt];  % Time from S to before R wave arrival
T3 = [tr];              % Arrival of R wave
T4 = [tr+dt:dt:tmax];   % Time after passage of R wave
T = [T0, T1, T2, T3, T4]; % Total time array
T = T*r/Cs;             % actual (physical) time

T12 = T1.^2;
T22 = T2.^2;
S11 = sqrt(T12-x1);
S21 = sqrt(T12-x2);
S31 = sqrt(x3-T12);
S12 = sqrt(T22-x1);
S22 = sqrt(T22-x2);
S32 = sqrt(x3-T22);

% I.- VERTICAL LOAD
% Vertical displacements due to vertical step load
f = (1-pois)/(2*pi*mu*r);
if (imag(x1)==0)
   A1 = (x1-0.5)^2*sqrt(a2-x1)/((x1-x2)*(x1-x3));
   A2 = (x2-0.5)^2*sqrt(a2-x2)/((x2-x1)*(x2-x3));
   A3 = (x3-0.5)^2*sqrt(x3-a2)/((x3-x1)*(x3-x2));
   U1 = 0.5*f*(1-A1./S11-A2./S21-A3./S31);
else
   A1 = (x1-0.5)^2*sqrt(a2-x1)/((x1-x2)*(x1-x3));
   A3 = (x3-0.5)^2*sqrt(x3-a2)/real((x3-x1)*(x3-x2));
   U1 = 0.5*f*(1-2*real(A1./S11)-A3./S31);
end
U2 = f*(1-A3./S32);
U3 = [-sign(A3)*inf];
U4 = f*ones(size(T4));
U4(1) = U2(length(U2));
Uzz = [U0, U1, U2, U3, U4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Radial displacements due to vertical load
n1 = b2/(a2-x1);
if pois==0
   x2 = a2;
   n2 = inf;
else
   n2 = b2/(a2-x2);
end
n3 = b2/(a2-x3);
k2 = (T12-a2)/b2;	% = k^2
if imag(x1)==0
  B1 = ellipint3(90,n1*k2,k2)*(1-2*x1)*(1-x1)/(x1-x2)/(x1-x3);
  B2 = ellipint3(90,n2*k2,k2)*(1-2*x2)*(1-x2)/(x2-x1)/(x2-x3);
  B3 = ellipint3(90,n3*k2,k2)*(1-2*x3)*(1-x3)/(x3-x1)/(x3-x2);
  U1 = 2*ellipint3(90,0,k2)-B1-B2-B3;
else
  B1 = 2*real(ellipint3(90,n1*k2,k2)*(1-2*x1)*(1-x1)/(x1-x2)/(x1-x3));
  tmp = real((x3-x1)*(x3-x2));
  B3 = ellipint3(90,n3*k2,k2)*(1-2*x3)*(1-x3)/tmp;
  U1 = 2*ellipint3(90,0,k2)-B1-B3;
end
fac = 1/(8*pi^2*mu*r);
f = fac/sqrt(b2^3);
U1 = f*U1.*T1;
t2 = T2.^2;
k2 = b2./(t2-a2);	% inverse of k^2
if imag(x1)==0
  B1 = ellipint3(90,n1,k2)*(1-2*x1)*(1-x1)/(x1-x2)/(x1-x3);
  B2 = ellipint3(90,n2,k2)*(1-2*x2)*(1-x2)/(x2-x1)/(x2-x3);
  B3 = ellipint3(90,n3,k2)*(1-2*x3)*(1-x3)/(x3-x1)/(x3-x2);
  U2 = 2*ellipint3(90,0,k2) - B1 - B2 - B3;
else
  B1 = 2*real(ellipint3(90,n1,k2)*(1-2*x1)*(1-x1)/(x1-x2)/(x1-x3));
  tmp = real((x3-x1)*(x3-x2));
  B3 = ellipint3(90,n3,k2)*(1-2*x3)*(1-x3)/tmp;
  U2 = 2*ellipint3(90,0,k2) - B1 - B3;
end
C = (2*x3-1)^3/(1-4*x3+8*b2*x3^3);
U2 = f*U2.*T2.*sqrt(k2);
U3 = U2(length(U2));
U4 = 2*pi*fac*C*T4./sqrt(T4.^2-x3);
Urz = [U0, U1, U2, U3, U4];

% II.- HORIZONTAL LOAD
% Radial displacements due to horizontal load
f = 1/(2*pi*mu*r);
fac = (1-pois)*f;
C1 = (1-x1)*sqrt(a2-x1)/((x1-x2)*(x1-x3));
if imag(x1)==0
  C2 = (1-x2)*sqrt(a2-x2)/((x2-x1)*(x2-x3));
  C3 = (1-x3)*sqrt(x3-a2)/((x3-x1)*(x3-x2));
  U1 = 0.5*fac*T12.*(C1./S11+C2./S21+C3./S31);
else
  C3 = (1-x3)*sqrt(x3-a2)/real((x3-x1)*(x3-x2));
  U1 = 0.5*fac*T12.*(2*real(C1./S11)+C3./S31);
end
U2 = f+fac*C3*T22./S32;
U3 = f;
U4 = f*ones(size(T4));
Urx = [U0, U1, U2, U3, U4];

% Tangential displacements due to horizontal load
if imag(x1)==0
  U1 = 0.5*fac*(1-C1*S11-C2*S21+C3*S31);
else
  U1 = 0.5*fac*(1-2*real(C1*S11)+C3*S31);
end
U2 = fac*(1+C3*S32);
U3 = fac; 
U4 = fac*ones(size(T4));
Utx = [U0, U1, U2, U3, U4];

% Displacements on epicentral axis
% ********************************
z = 1;  % depth at which displacements are computed
t = [T2,T3,T4]; % time after arrival of S wave
t2 = t.^2;
D1 = (16*(1-a2)*t2+8*(3-8*a2+6*a2^2)).*t2;
D1 = (D1+8*(1-6*a2+10*a2^2-6*a2^3)).*t2+(1-2*a2)^4;
D2 = (16*(1-a2)*t2-8*(3-4*a2)).*t2;
D2 = (D2+8*(1-2*a2)).*t2+1;
S1 = 0.5*(1+sqrt(1+(1-a2)./t2)).*D1;
S2 = 0.5*(1+sqrt(1-(1-a2)./t2)).*D2;

% a) Horizontal load
fac = 0.25/pi/mu/z;
U1 = sqrt(T12-a2+1);
U1 = 2*T1.*(T12-a2).*U1;
U1 = U1./(U1-(2*T12-2*a2+1).^2);
U2 = ((128*(1-a2)*t2-64*(1+4*a2-6*a2^2)).*t2-16*(3-15*a2-4*a2^2+24*a2^3)).*t2;
U2 = ((U2+16*a2*(4-17*a2+10*a2^2+8*a2^3)).*t2+16*a2*(1-3*a2+7*a2^2-6*a2^3)).*t2;
U2 = (U2-(1-10*a2+40*a2^2-48*a2^3+16*a2^4)).*t2./D1./D2+1;
U2 = U2-(1-a2)*((t2-a2).*(2*t2-2*a2+1).^2./S1+2*t2.*(2*t2-1).*(t2-1)./S2);
Vxx = fac*[U0,U1,U2];   % Horiz. displacement at depth due to horiz. load

% b) Vertical load
fac = 0.5/pi/mu/z;
U1 = 2*T12-2*a2+1;
U1 = T12.*U1./(U1.^2-4*T1.*(T12-a2).*sqrt(T12-a2+1));
U2 = ((128*(1-a2)*t2+64*(1-2*a2)*(2-4*a2+a2^2)).*t2-16*(21-37*a2+4*a2^2+36*a2^3-16*a2^4)).*t2;
U2 = ((U2+16*(3+26*a2-78*a2^2+70*a2^3-8*a2^4-8*a2^5)).*t2+4*(15-87*a2+116*a2^2+24*a2^3-136*a2^4+64*a2^5)).*t2;
U2 = (U2+(11-28*a2+16*a2^2)*(1-2*a2)^3).*t2./D1./D2;
U2 = U2+(1-a2)*( 2*t2.*(t2-a2).*(2*t2-2*a2+1)./S1 + (2*t2-1).^2.*(t2-1)./S2);
Vzz = fac*[U0,U1,U2];    % Vert. displacement at depth due to horiz. load

if NoPlot
  % omit all plots
  if nargout==0, 
   clear T Urx Utx Uzz Urz
  end
  return
end

% Plot response functions
% ***********************

plot(T, Uzz);
grid on;
axis ([0 tmax -1 0.4]);
tit = sprintf( ... 
   'Vertical displacement due to vertical point (step) load, \\nu=%5.3f',pois);
title (tit);
titx = 'Dimensionless time';
xlabel(titx);
EasyPlot('Uzz', tit, titx, T, Uzz, 'r', pois);
pause;

plot (T,Urz);
grid on;
axis ([0 tmax -0.2 0.6]);
tit = sprintf( ...
  'Radial displacement due to vertical point (step) load, \\nu=%5.3f',pois);
title (tit);
xlabel(titx);
EasyPlot('Urz', tit, titx, T, Urz, 'r', pois);
pause;

plot (T,Urx);
grid on;
tit = sprintf ( ...
  'Radial displacement due to horizontal point (step) load, \\nu=%5.3f',pois);
title (tit);
xlabel(titx);
EasyPlot('Urx', tit, titx, T, Urx, 'r', pois);
pause;

plot (T,Utx);
grid on;
tit = sprintf( ...
  'Tangential displacement due to horizontal point (step) load, \\nu=%5.3f',pois);
title (tit);
xlabel(titx);
EasyPlot('Utx', tit, titx, T, Utx, 'r', pois);
pause;

plot (T,Vzz);
grid on;
tit = sprintf( ...
  'Vertical displacement under load at epicentral axis, \\nu=%5.3f',pois);
title(tit);
xlabel(titx);
EasyPlot('Axis_zz', tit, titx, T, Vzz, 'r', pois);
pause;

plot (T,Vxx);
grid on;
tit = sprintf( ...
  'Horizontal displacement under load at epicentral axis, \\nu=%5.3f',pois);
title(tit);
xlabel(titx);
EasyPlot('Axis_xx', tit, titx, T, Vxx, 'r', pois);
pause;

close all

if nargout==0, 
 clear T Urx Utx Uzz Urz
end
return


function EasyPlot(fname, tit, titx, T, U, flag, pois)
% Creates files for display with the EasyPlot program
nu = floor(1000*pois);
if nu>99
  fname = sprintf('%s%3d.ezp',fname,nu);
elseif nu>9
  fname = sprintf('%s0%2d.ezp',fname,nu);
else
  fname = sprintf('%s00%1d.ezp',fname,nu);
end
fout = fopen (fname, 'w');
fprintf (fout, '/et g "%s"\n', tit);
fprintf (fout, '/et x "%s"\n', titx);
fprintf (fout, '/og on\n');
fprintf (fout, '/sd off\n');
fprintf (fout, '/sm off\n');
[nr,nc] = size(U);
for n=1:nr
   fprintf (fout, '//nc\n');
   if flag=='r'
	  fprintf (fout, '%15.5f  %15.5f\n', [T;U(n,:)]);
   elseif flag=='c'
	  fprintf (fout, '%15.5f  %15.5f  %15.5f\n', [T;real(U(n,:));imag(U(n,:))]);
   else
   'Not a valid option in EasyPlot';
   end
end 
fclose (fout);
return


function [el3]=ellipint3(phi,N,M);
%    [EL3] = ELLIPINT3 (phi,N,M) returns the elliptic integral of
%            the third kind, evaluated for each value of N, M
%            Can also be used to obtain the elliptic integral
%            of the first kind by setting N=0.
%    Arguments :  phi   -- Upper limit of integration (in degrees)
%                 M=[m] -- Modulus   (some authors use k=sqrt(m))
%                          M can be a scalar or vector
%                 N=[n] -- Parameter (some authors use c=-n)
%                          N can be also be scalar or vector, but if
%                          the latter, it must agree in size with M
%    Definition: If n, m are elements of N, M, then
%
%                phi
%    el3 = integ  |   [dt/((1+n*sin^2(t))*sqrt(1-m*sin^2(t)))
%                 0
%
%    Observe that m = k^2 is the square of the argument
%    used by some authors for the elliptic integrals
%    Method:  10-point Gauss-Legendre quadrature
if (phi < 0)
  'Error, first argument in ellipint3 cannot be negative'
  el3 = 0;
  return
end
if length(N) ==1
  N = N*ones(size(M));
elseif length(N) ~= length(M)
  'Error, wrong size of second argument in ellipint3.'
  'Should be size(N)=1 or size(N)=size(M)'
  el3 = 0;
  return
end
tol = 1.e-8;
ang = phi*pi/180;
psi = ang/2;
t = [.9931285991850949,.9639719272779138,.9122344282513259, ...
  .8391169718222188,.7463319064601508,.6360536807265150, ...
  .5108670019508271,.3737060887154195,.2277858511416451, ...
  .7652652113349734d-1];
w = [.1761400713915212d-1,.4060142980038694d-1,.6267204833410907d-1, ...
  .8327674157670475d-1,.1019301198172404,.1181945319615184, ...
  .1316886384491766,.1420961093183820,.1491729864726037,.1527533871307258];
t1 = psi*(1+t);
t2 = psi*(1-t);
s1 = sin(t1).^2;
s2 = sin(t2).^2;
el3 = zeros(size(M));
s = sin(ang)^2;
for j=1:length(M)
  k2 = M(j);
  n = N(j);
  % assuming phi is in degrees here
  if (phi<=90 & abs(1+n*s)>tol & abs(1-k2*s)>tol) | ...
      (phi>90 & abs(1+n)>tol & abs(1-k2)>tol)  
    f1 = 1./( (1+n*s1).*sqrt(1-k2*s1) );
    f2 = 1./( (1+n*s2).*sqrt(1-k2*s2) );
    el3(j) = (f1+f2)*w';
  else
    el3(j) = inf;
  end
end
el3 = psi*el3;
return
