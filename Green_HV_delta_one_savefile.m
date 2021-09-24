
function Ky=Green_HV_delta_one_savefile(Cs,Cp,Rho,H,Damp_s,Damp_p,Ntps,dt,ymax,delay_t,fmax2,flag,infold);

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
%
%  Ntps   = total number of time step
%  dt = time step
%  y_pos = receiver distance (row), 0 is forbidden
%  delay_t   = excitation time of delta function δ(t-delay_t)
%  flag  = 'V' for vertical load and 'H' for horizontal load
%  fmax2 = max. frequency to be calculated
% In addition, the user must specify a material properties of the soil,
% which are defined in the block of statements immediately below.
% Material properties (size of arrays = number of layers)
%**************************************************************************
% Cs   = [300 180 450];   % Shear wave velocity [m/s]
% Cp   =[600 1700 2000];
% Rho  = [1 1 1 ]*2000;    % Mass density of layers (kg/m^3)
% H    = [20 20 inf];      % depth of layers [m]
% Damp_s=1/(2*Qs);  [0.02 0.02 0.02]  %  viscoealstic damping
% Damp_p=1/(2*Qp);  [0.02 0.02 0.02]  % A zero represents elastic media.

% Cs   = [200 400];   % Shear wave velocity [m/s]
% Cp   =[800 1200];
% Rho  = [1 1 ]*2000;    % Mass density of layers (kg/m^3)
% H    = [10 inf];      % depth of layers [m]

% If the depth of the last layer is infinite, then it is a half-space,
% otherwise it is a stratum (layer over fixed base)
% Related material properties
% Cp = Cs.*sqrt((2-2*Pois)/(1-2*Pois)); % P-wave velocity
%**************************************************************************
% output arguments
% Wavef(:,:,1)-Radial compenent
% Wavef(:,:,2)-Tangential compenent
% Wavef(:,:,3)-Vertical compenent
%**************************************************************************
% example: homogenous half space
% tic;[ut3,t3]=Kausel_Kai_Green_HV([400],[1200],2,[inf],[0],[0],4000,0.001,50:50:1000,10,'V');toc
%**************************************************************************

% Ntps=2^nextpow2(Ntps);  %% There is no need for fft and ifft in matlab with 2^n
%% Default input parameters
%ymax=max(y_pos);   % max. range of Y coordinate
if mod(Ntps,2)>0
    Ntps=Ntps-1;
end
Nperiod=ceil(Ntps*dt*max(Cp)/ymax)+1; % make sure t<sqrt((Ly-ymax)^2+z^2)/max(Cp)
Ly = Nperiod*ymax;      % spatial length for periodicity in Y
t=0:dt:(Ntps-1)*dt; % time sequence

fmin=0;
fmax=1/dt/2;
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
kmax=sqrt(2)*max(kmax1,kmax2);
%kmax=w(end)/min(Cs);
%**************************************************************************
% The following step make a long trouble for me since it causes a wavelet
% in the just exciting time for every traces.
%**************************************************************************

dky = 2*pi/Ly;  % wavenumber step in Y
nky = ceil(kmax/dky);

if mod(nky,2)>0, nky=nky+1; end % make nky even
Ky = 0:dky:nky*dky; % wavenumber vector

% Ricker wavelet
delta_w=2*pi*df;  %% or 3*pi*df

% source_tt=rickerkai(t,fm).*exp(-delta_w*t); %%  Exponential Window Method
% arrf=(fft(source_tt));

%infold='D:\Modeling noise\';
if ~exist(infold,'dir')
    mkdir(infold)
end
delete(strcat(infold,'*.mat'));

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
    
    w_ew=w(ii)-1i*delta_w; %% Exponential Window Method  compare with w_ew=w(ii);
    UA1=zeros(length(Ky),1);
    UA2=zeros(length(Ky),1);
    UA3=zeros(length(Ky),1);
           
    if ~any(strcmp(flag,{'H','h','V','v'}))
        error('Wrong input strings within ({H,h,V,v})');
    end
    %**************************************************************************
    if any(strcmp(flag,{'H','h'}))
        [KGsh, KGsvp]=arrayfun(@(ky)(GlobalStiffMat(Cp,Cs,Damp_p,Damp_s,Rho,H,ky,w_ew)),Ky,'UniformOutput',false);
        UA_temp_z=cellfun(@(x) x\qR,KGsvp,'UniformOutput',false);
        %   UA_temp=cellfun(@(x) pinv(x)*qR,KGsvp,'UniformOutput',false); %% pinv is not recommended,time consuming
        UA_z=cell2mat(UA_temp_z);
        
        UA_temp_t=cellfun(@(x) x\qT,KGsh,'UniformOutput',false);
        UA_t=cell2mat(UA_temp_t);
        
        if ~any(isnan(UA_z))
            UA1 = UA_z(1,:);
            UA1=UA1.';
            UA3 = UA_z(2,:);
            UA3=UA3.';
        end
        if ~any(isnan(UA_t))
            UA2 = UA_t(1,:);
            UA2=UA2.';
        end     
    end
    %**************************************************************************
    if any(strcmp(flag,{'V','v'}))
        KGsvp=arrayfun(@(ky)(GlobalStiffMat_R(Cp,Cs,Damp_p,Damp_s,Rho,H,ky,w_ew)),Ky,'UniformOutput',false);
        UA_temp_z=cellfun(@(x) x\qV,KGsvp,'UniformOutput',false);
        %   UA_temp=cellfun(@(x) pinv(x)*qR,KGsvp,'UniformOutput',false);%% pinv is not recommended,time consuming
        UA_z=cell2mat(UA_temp_z);
        
        if ~any(isnan(UA_z))
            UA1 = UA_z(1,:);
            UA1=UA1.';
            UA3 = UA_z(2,:);
            UA3=UA3.';
        end
    end
    ukai=[UA1 UA2 UA3];
    filename=strcat(num2str(ii),'.mat');
    parsave([infold,filename],ukai);

    pause(pauseTime);
    % increment counter to track progress
    ppm.increment();
end
toc
% Delete the progress handle when the parfor loop is done. 
delete(ppm);

