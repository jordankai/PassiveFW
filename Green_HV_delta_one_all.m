
function UK=Green_HV_delta_one_all(Cs,Cp,Rho,H,Damp_s,Damp_p,Ntps,dt,ymax,delay_t,fmax2);

% The code is based on E.Kausel's moving bell load code
%**************************************************************************
%    (many thanks to Professor E.Kausel)
%    Written and Copyleft by Eduardo Kausel, MIT, kausel@mit.edu
%    Version 1, May 21, 2017
%    You may use and modify this code, provided you acknowledge the source
% Referrence "Generalized stiffness matrix method for layered soils" (2018)
%**************************************************************************
%
% Written by Zhang Kai for calculating 3D green's function in time domain
% using discrete wavenumber method.
% 2021/6/26  naturekai@126.com
%
%**************************************************************************
%
% Computes the response(Green's function) in space and time due to a narrowly
% peaked 3-D bell load acting at the upper free surface of a layered elastic
% and viscoelastic half-space.
% The output is given at the surface only
%
% Input arguments
%**************************************************************************
%  Ntps   = total number of time step
%  dt = time step
%  ymax = max receiver distance, 0 is forbidden
%  delay_t   = excitation time of delta function δ(t-delay_t)
%  fmax2 = max. frequency to be calculated
% In addition, the user must specify a material properties of the soil,
% which are defined in the block of statements immediately below.
% Material properties (size of arrays = number of layers)
%**************************************************************************
% If the depth of the last layer is infinite, then it is a half-space,
% otherwise it is a stratum (layer over fixed base)
% Related material properties
% Cp = Cs.*sqrt((2-2*Pois)/(1-2*Pois)); % P-wave velocity
%**************************************************************************
% output arguments
% UK.Ky- Wavenumber vector
% UK.ukai1-Radial compenent due to Horizontal point load
% UK.ukai2-Tangential compenent due to Horizontal point load
% UK.ukai3-Vertical compenent due to Horizontal point load
% UK.ukai4-Radial compenent due to vertical point load
% UK.ukai5-Vertical compenent due to vertical point load
%**************************************************************************

% Ntps=2^nextpow2(Ntps);  %% There is no need for fft and ifft in matlab with 2^n
if mod(Ntps,2)>0
    Ntps=Ntps-1;
end
Nperiod=ceil(Ntps*dt*max(Cp)/ymax)+1; % make sure t<sqrt((Ly-ymax)^2+z^2)/max(Cp)
Ly = Nperiod*ymax;      % spatial length for periodicity in Y
fmin=0; % min frequency
fmax=1/dt/2; % max frequency
if fmax2>fmax
    fmax2=fmax;
end
ff=linspace(fmin,fmax,Ntps/2+1);
df=ff(2)-ff(1);% frequency step
fN=ceil(fmax2/df)+1;
w=2*pi*ff;

nL = length(Cs);             % Number of layers (materials)
nL2 = 2*nL;                  % double the number of layers

kmax1=2*w(end)/min(Cs);  % max wavenumber according to the article Kausel(2018)
kmax2=sqrt((1.5*pi/min(H))^2+(w(end)/min(Cs))^2);

% kmax1=2*2*pi*fmax2/min(Cs);  % max wavenumber according to the article Kausel(2018)
% kmax2=sqrt((1.5*pi/min(H))^2+(2*pi*fmax2/min(Cs))^2);

kmax=sqrt(2)*max(kmax1,kmax2);
dky = 2*pi/Ly;  % wavenumber step in Y
nky = ceil(kmax/dky);
if mod(nky,2)>0, nky=nky+1; end % make nky even
Ky = 0:dky:nky*dky; % wavenumber vector

% Exponential Window Method
delta_w=2*pi*df;  %% or 3*pi*df

ukai1=zeros(nky+1,fN);
ukai2=zeros(nky+1,fN);
ukai3=zeros(nky+1,fN);
ukai4=zeros(nky+1,fN);
ukai5=zeros(nky+1,fN);

numIterations = fN;
cpu_N=feature('numCores');
ppm = ParforProgressbar(numIterations,'parpool', {'local', cpu_N});
pauseTime = 10/numIterations;
tic
parfor ii=2:fN %% loops over the frequency
    % do fft integration using delta function with f(t): = f(0)
    % f=exp(-delta_w*t).* exp(-1i*w*t), f(0)=1
    factork=exp(-1i*w(ii)*delay_t)*exp(-delta_w*delay_t);  %% Time delay(0.5s) so that separate green with wavelet
    qR = zeros(nL2,1);  qR(1)=factork/2/pi;        % load vector for Radial
    qT = zeros(nL,1);   qT(1)=factork/2/pi;        % load vector for Tangential
    qV = zeros(nL2,1);  qV(2)=factork/2/pi;        % load vector for Vertical
    
    w_ew=w(ii)-1i*delta_w; %% Exponential Window Method
    UA1=zeros(length(Ky),1);
    UA2=zeros(length(Ky),1);
    UA3=zeros(length(Ky),1);
    UA4=zeros(length(Ky),1);
    UA5=zeros(length(Ky),1);
    %**************************************************************************
    
    [KGsh, KGsvp]=arrayfun(@(ky)(GlobalStiffMat(Cp,Cs,Damp_p,Damp_s,Rho,H,ky,w_ew)),Ky,'UniformOutput',false);
    UA_temp=cellfun(@(x) x\qR,KGsvp,'UniformOutput',false);
    %   UA_temp=cellfun(@(x) pinv(x)*qR,KGsvp,'UniformOutput',false); %% pinv is not recommended,time consuming
    UA=cell2mat(UA_temp);            
    if ~any(isnan(UA))
        UA1 = UA(1,:);
        UA1=UA1.';
        UA3 = UA(2,:);
        UA3=UA3.';
    end
    
    UA_temp=cellfun(@(x) x\qT,KGsh,'UniformOutput',false);
    UA=cell2mat(UA_temp);
    if ~any(isnan(UA))
        UA2 = UA(1,:);
        UA2=UA2.';
    end   
    %**************************************************************************
    
    UA_temp=cellfun(@(x) x\qV,KGsvp,'UniformOutput',false);
    %   UA_temp=cellfun(@(x) pinv(x)*qR,KGsvp,'UniformOutput',false);%% pinv is not recommended,time consuming
    UA=cell2mat(UA_temp);    
    if ~any(isnan(UA))
        UA4 = UA(1,:);
        UA4=UA4.';
        UA5 = UA(2,:);
        UA5=UA5.';
    end
    
    ukai1(:,ii)=UA1;
    ukai2(:,ii)=UA2;
    ukai3(:,ii)=UA3;
    ukai4(:,ii)=UA4;
    ukai5(:,ii)=UA5;
    
    pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end
toc

% Delete the progress handle when the parfor loop is done.
delete(ppm);
UK.k=Ky;
UK.HR=ukai1;
UK.HT=ukai2;
UK.HV=ukai3;
UK.VR=ukai4;
UK.VV=ukai5;
