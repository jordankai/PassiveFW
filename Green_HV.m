
function [Wavef,t]=Green_HV(Cs,Cp,Rho,H,Damp_s,Damp_p,Ntps,dt,y_pos,fm,flag);

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
%  fm   = main frequency of source wavelet
%  flag  = 'V' for vertical load and 'H' for horizontal load

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

% Default input parameters
ymax=max(y_pos);   % max. range of Y coordinate
if mod(Ntps,2)>0
    Ntps=Ntps-1;
end
Nperiod=ceil(Ntps*dt*max(Cp)/ymax)+1; % make sure t<sqrt((Ly-ymax)^2+z^2)/max(Cp)
Ly = Nperiod*ymax;      % spatial length for periodicity in Y
N_max=floor(((Ly-ymax)/max(Cp))/dt);
if Ntps>N_max
   Ntps=N_max;
end

if mod(Ntps,2)>0
    Ntps=Ntps-1;
end
t=0:dt:(Ntps-1)*dt; % time sequence

fmin=0;
fmax=1/dt/2;
ff=linspace(fmin,fmax,Ntps/2+1);
df=ff(2)-ff(1);% frequency step
w=2*pi*ff;

nL = length(Cs);             % Number of layers (materials)
nL2 = 2*nL;                  % double the number of layers

kmax1=2*w(end)/min(Cs);  % max wavenumber according to the article Kausel(2018)
kmax2=sqrt((1.5*pi/min(H))^2+(w(end)/min(Cs))^2);
kmax=4*max(kmax1,kmax2);
%kmax=w(end)/min(Cs);
%**************************************************************************
% The following step make a long trouble for me since it causes a wavelet
% in the just exciting time for every traces.
%**************************************************************************
Y_num=length(y_pos);
dky = 2*pi/Ly;  % wavenumber step in Y
nky = ceil(kmax/dky);
if nky<Y_num
    nky = Y_num;
end
if mod(nky,2)>0, nky=nky+1; end % make nky even
Ky = 0:dky:nky*dky; % wavenumber vector

% Ricker wavelet
delta_w=2*pi*df;  %% or 2*pi*df

source_tt=rickerkai(t,fm).*exp(-delta_w*t); %%  Exponential Window Method
arrf=(fft(source_tt));

% Note: The reason for using a bell load instead of a point load is that
% the bell decays rather fast in the wavenumber domain, yet it has
% virtually the same effect as the point load.

r0 = 0.02*y_pos;  %% make sure bell load has approximately the same effect as the point load
r0(r0>1)=1;
r2 = r0.^2;
%Ex = exp(-(Kx2*r2/4));       % ~Hankel transform of Gaussian bell
Ey = exp(-(Ky.^2.*r2'/4)); % ~Hankel transform of Gaussian bell

%**************************************************************************
%  green's function initialization
ukai1=zeros(Y_num,length(w));
ukai2=zeros(Y_num,length(w));
ukai3=zeros(Y_num,length(w));
ukai=zeros(Y_num,length(w),3);
Wavef=zeros(Ntps,Y_num,3);
% Bessel function for order 0 and 1
besl_kai0=besselj(0,Ky'.*y_pos);
besl_kai1=besselj(1,Ky'.*y_pos);

ParN = (length(w)-1);  %% parallel progress monitor
parfor_progress(ParN);
parfor ii=2:(length(w)-1) %% loops over the frequency
    
    qR = zeros(nL2,1);  qR(1)=arrf(ii);        % load vector for Radial
    qT = zeros(nL,1);   qT(1)=arrf(ii);        % load vector for Tangential
    qV = zeros(nL2,1);  qV(2)=arrf(ii);        % load vector for Vertical
    
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
            UA1 = UA_z(1,:).*Ey;
            UA1=UA1.';
            UA3 = UA_z(2,:).*Ey;
            UA3=UA3.';
        end
        if ~any(isnan(UA_t))
            UA2 = UA_t(1,:).*Ey;
            UA2=UA2.';
        end
        
        % inverse Hankel transform in Y via plain summation
        UA_sum1= (Ky'.*UA1.*besl_kai0-UA1.*besl_kai1./y_pos+UA2.*besl_kai1./y_pos);
        ukai1(:,ii)=(sum(UA_sum1(2:end-1,:))+0.5*(UA_sum1(end,:)+UA_sum1(1,:)))*dky; %% trapezoid rule(Discrete wavenumber method)
        
        UA_sum2= (Ky'.*UA2.*besl_kai0-UA2.*besl_kai1./y_pos+UA1.*besl_kai1./y_pos);
        ukai2(:,ii)=(sum(UA_sum2(2:end-1,:))+0.5*(UA_sum2(end,:)+UA_sum2(1,:)))*dky; %% trapezoid rule
        
        UA_sum3= Ky'.*UA3.*besl_kai1;
        ukai3(:,ii)=(sum(UA_sum3(2:end-1,:))+0.5*(UA_sum3(end,:)+UA_sum3(1,:)))*dky; %% trapezoid rule
    end
    %**************************************************************************
    if any(strcmp(flag,{'V','v'}))
        KGsvp=arrayfun(@(ky)(GlobalStiffMat_R(Cp,Cs,Damp_p,Damp_s,Rho,H,ky,w_ew)),Ky,'UniformOutput',false);
        UA_temp_z=cellfun(@(x) x\qV,KGsvp,'UniformOutput',false);
        %   UA_temp=cellfun(@(x) pinv(x)*qR,KGsvp,'UniformOutput',false);%% pinv is not recommended,time consuming
        UA_z=cell2mat(UA_temp_z);
        
        if ~any(isnan(UA_z))
            UA1 = UA_z(1,:).*Ey;
            UA1=UA1.';
            UA3 = UA_z(2,:).*Ey;
            UA3=UA3.';
        end
        % inverse Hankel transform in Y via plain summation
        UA_sum1= -Ky'.*UA1.*besl_kai1;
        ukai1(:,ii)=(sum(UA_sum1(2:end-1,:))+0.5*(UA_sum1(end,:)+UA_sum1(1,:)))*dky; %% trapezoid rule
        
        UA_sum3= Ky'.*UA3.*besl_kai0;
        ukai3(:,ii)=(sum(UA_sum3(2:end-1,:))+0.5*(UA_sum3(end,:)+UA_sum3(1,:)))*dky; %% trapezoid rule
    end
    parfor_progress;
end
parfor_progress(0);
% ukai(:,:,1)=ukai1*cos(az_ang);
% ukai(:,:,2)=-ukai2*sin(az_ang);
% ukai(:,:,3)=ukai3*cos(az_ang);
ukai(:,:,1)=ukai1;
ukai(:,:,2)=ukai2;
ukai(:,:,3)=ukai3;
for ii=1:3
%     trace_f=ukai(:,:,ii); 
    trace_f=Haning_kai(ukai(:,:,ii)); %% Haning window
    % exp(1i*pi/4) Adjusting the phase difference (pi/4)
    % between 2D and 3D green's functions to ensure that they are consistent
    % trace_f=costap_filter(ukai(:,:,ii),0.1,2); %%%% cosine windows in frequency domain is not suitable
    % inverse fourier transform to time domain
    trace_t=[trace_f,  conj(trace_f(:,end-1:-1:2))];
    ut=(ifft(trace_t,[],2));
    ut=ut.*exp(delta_w*t);  %% Exponential Window Method
    ut=ut.';
  %  ut=mute_kai(ut,1.5*max(Cp),y_pos,dt); %% mute events faster than 1.5*max(Cp)
    Wavef(:,:,ii)=ut;
end
end
% fmat=moviein(512);
% for j=70:512;
% surf(Wave(8:8:end,2:2:end,j)/max(max(abs(Wave(8:8:end,2:2:end,j)))))
% view(80,20)
% fmat(:,j)=getframe;
% end
% movie(fmat,1)
