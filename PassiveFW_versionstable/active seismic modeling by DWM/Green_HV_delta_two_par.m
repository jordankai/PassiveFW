
function [Wavef,t]=Green_HV_delta_two_par(UK,Ntps,dt,y_pos,fmax2,flag);

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
% 2021/9/30  naturekai@126.com
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
%  UK- structure data in Wavenumber domain
%  Ntps   = total number of time step
%  dt = time step
%  y_pos = receiver distance (row), 0 is forbidden
%  fmax2 = max. frequency to be calculated
%  flag  = 'V' for vertical load and 'H' for horizontal load
%          'all' or 'All' for both of them

%**************************************************************************
% output arguments
% Wavef(:,:,1)-Radial compenent due to Horizontal point load
% Wavef(:,:,2)-Tangential compenent due to Horizontal point load
% Wavef(:,:,3)-Vertical compenent due to Horizontal point load
% Wavef(:,:,4)-Tangential compenent due to vertical point load
% Wavef(:,:,5)-Vertical compenent due to vertical point load
% t-time sequence
%**************************************************************************

% Ntps=2^nextpow2(Ntps);  %% There is no need for fft and ifft in matlab with 2^n
%% Default input parameters
if ~any(strcmp(flag,{'H','h','V','v','hv','HV','Hv','hV'}))
    error('Wrong input strings within ({H,h,V,v,hv,HV,Hv,hV})');
end

Ky=UK.k;
uk1=UK.HR;
uk2=UK.HT;
uk3=UK.HV;
uk4=UK.VR;
uk5=UK.VV;

if mod(Ntps,2)>0
    Ntps=Ntps-1;
end
t=0:dt:(Ntps-1)*dt; % time sequence
dky=Ky(2)-Ky(1);
fmin=0;
fmax=1/dt/2;
if fmax2>fmax
    fmax2=fmax;
end
ff=linspace(fmin,fmax,Ntps/2+1);
df=ff(2)-ff(1);% frequency step
w=2*pi*ff;
fN=ceil(fmax2/df)+1;
Y_num=length(y_pos);

% Exponential Window Method
delta_w=2*pi*df;  %% or 3*pi*df

% Note: The reason for using a bell load instead of a point load is that
% the bell decays rather fast in the wavenumber domain, yet it has
% virtually the same effect as the point load.

r0 = 0.02*y_pos;  %% make sure bell load has approximately the same effect as the point load
r0(r0>1)=1;
r2 = r0.^2;
%Ex = exp(-(Kx2*r2/4));       % ~Hankel transform of Gaussian bell
Ey = exp(-(Ky'.^2.*r2/4)); % ~Hankel transform of Gaussian bell

%**************************************************************************
%  green's function initialization
ukai1=zeros(Y_num,length(w));
ukai2=zeros(Y_num,length(w));
ukai3=zeros(Y_num,length(w));
ukai4=zeros(Y_num,length(w));
ukai5=zeros(Y_num,length(w));
ukai=zeros(Y_num,length(w),3);
Wavef=zeros(Ntps,Y_num,3);
% Bessel function for order 0 and 1
besl_kai0=besselj(0,Ky'.*y_pos);
besl_kai1=besselj(1,Ky'.*y_pos);

numIterations = fN;
cpu_N=feature('numCores');
str_kai=strcat(['Parallel computation using ', num2str(cpu_N)], ' cpu cores, Please wait...');
hbar = parfor_progressbar(numIterations,str_kai); %create the progress bar

parfor ii=2:fN %% loops over the frequency
    
    %**************************************************************************
    if any(strcmp(flag,{'H','h','hv','HV','Hv','hV'}))
        
        % inverse Hankel transform in Y via plain summation
        UA1=uk1(:,ii).*Ey;
        UA2=uk2(:,ii).*Ey;
        UA_sum1= (Ky'.*UA1.*besl_kai0-UA1.*besl_kai1./y_pos+UA2.*besl_kai1./y_pos);
        ukai1(:,ii)=(sum(UA_sum1(2:end,:))+0.5*UA_sum1(1,:))*dky; %% (Discrete wavenumber method)
        % ukai1(:,ii)=(sum(UA_sum1(2:end-1,:))+0.5*(UA_sum1(end,:)+UA_sum1(1,:)))*dky; %% trapezoid rule
        
        UA_sum1= (Ky'.*UA2.*besl_kai0-UA2.*besl_kai1./y_pos+UA1.*besl_kai1./y_pos);
        ukai2(:,ii)=(sum(UA_sum1(2:end,:))+0.5*UA_sum1(1,:))*dky; %% (Discrete wavenumber method)
        % ukai2(:,ii)=(sum(UA_sum1(2:end-1,:))+0.5*(UA_sum1(end,:)+UA_sum1(1,:)))*dky; %% trapezoid rule
        
        UA1=uk3(:,ii).*Ey;
        UA_sum1= Ky'.*UA1.*besl_kai1;
        ukai3(:,ii)=(sum(UA_sum1(2:end,:))+0.5*UA_sum1(1,:))*dky; %% (Discrete wavenumber method)
        % ukai3(:,ii)=(sum(UA_sum1(2:end-1,:))+0.5*(UA_sum1(end,:)+UA_sum1(1,:)))*dky; %% trapezoid rule
        
    end
    %**************************************************************************
    if any(strcmp(flag,{'V','v','hv','HV','Hv','hV'}))
        
        % inverse Hankel transform in Y via plain summation
        UA1=uk4(:,ii).*Ey;
        UA_sum1= -Ky'.*UA1.*besl_kai1;
        ukai4(:,ii)=(sum(UA_sum1(2:end,:))+0.5*UA_sum1(1,:))*dky; %% (Discrete wavenumber method)
        % ukai4(:,ii)=(sum(UA_sum1(2:end-1,:))+0.5*(UA_sum1(end,:)+UA_sum1(1,:)))*dky; %% trapezoid rule
        
        UA1=uk5(:,ii).*Ey;
        UA_sum1= Ky'.*UA1.*besl_kai0;
        ukai5(:,ii)=(sum(UA_sum1(2:end,:))+0.5*UA_sum1(1,:))*dky; %% (Discrete wavenumber method)
        % ukai5(:,ii)=(sum(UA_sum1(2:end-1,:))+0.5*(UA_sum1(end,:)+UA_sum1(1,:)))*dky; %% trapezoid rule
    end
    
hbar.iterate(1); % update progress by one iteration
end
close(hbar); % close the progress bar

ukai(:,:,1)=ukai1;
ukai(:,:,2)=ukai2;
ukai(:,:,3)=ukai3;
ukai(:,:,4)=ukai4;
ukai(:,:,5)=ukai5;
for ii=1:5
    % trace_f=ukai(:,:,ii);
    trace_f=Haning_kai(ukai(:,:,ii)); %% Haning window
    % inverse fourier transform to time domain
    trace_t=[trace_f,  conj(trace_f(:,end-1:-1:2))];
    ut=(ifft(trace_t,[],2));
    ut=ut.*exp(delta_w*t);  %% Exponential Window Method
    ut=ut.';
    Wavef(:,:,ii)=ut;
end
end

