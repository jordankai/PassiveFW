%%  half-space model
clear;
tic
% Damp_p=1./(2*Qp); % damping
% Damp_s=1./(2*Qs);

model.vs   = 440;   % Shear wave velocity [m/s]
model.vp   =760;
model.dns  = 2;    % Mass density of layers (kg/m^3)
model.thk  = inf;      % depth of layers [m]
model.Damp_s=0;  %  viscoealstic damping
model.Damp_p=0;  % A zero represents elastic media.

%%  three-layer model
% model.vs   = [300 180 450];   % Shear wave velocity [m/s]
% model.vp   =[600 1700 2000];
% model.dns  = [2 2 2 ];    % Mass density of layers (kg/m^3)
% model.thk  = [20 20 inf];      % depth of layers [m]
% model.Damp_s=[0 0 0];  %  viscoealstic damping
% model.Damp_p=[0 0 0];  % A zero represents elastic media.

%% two-layer model
% model.vs   = [220 440];   % Shear wave velocity [m/s]
% model.vp   =[380 760];
% model.dns  = [1 1 ]*2000;    % Mass density of layers (kg/m^3)
% model.thk    = [10 inf];      % depth of layers [m]
% Qs=[30 30];
% Qp=[30 30];
% model.Damp_p =[0 0] ;%1./(2*Qp); % damping
% model.Damp_s = [0 0];%1./(2*Qs);

%% time and source setting
source.number=4; 
source.recordlength=10; % record time(s)
source.dt=0.001; %  time step
source.maxfre=100;  %% max frequency to be calculated
source.fmin=5;  %% min. Dominant Frequency of Ricker wavelet
source.fmax=30;  %% max. Dominant Frequency of Ricker wavelet

%%%%%%%%%%%%  set point loads   %%%%%%%%%%%%%%%%%%%%%
flag='h';

if any(strcmp(flag,{'v','V'}))
    source.Fx=zeros(1,source.number); % point load in x-axis
    source.Fy=zeros(1,source.number); % point load in y-axis
    source.Fz=-rand(1,source.number);% point load in z-axis,z-axis is positive in the upward direction
elseif any(strcmp(flag,{'h','H'}))
    source.Fx=zeros(1,source.number); % point load in x-axis
    source.Fy=rand(1,source.number); % point load in y-axis
    source.Fz=zeros(1,source.number);% point load in z-axis,z-axis is positive in the upward direction
elseif any(strcmp(flag,{'hv','HV','Hv','hV'}))
    source.Fx=rand(1,source.number); % point load in x-axis
    source.Fy=rand(1,source.number); % point load in y-axis
    source.Fz=-rand(1,source.number);% point load in z-axis,z-axis is positive in the upward direction
end

%% Geometry setting
SR=S_R_geometry(source.number);

%% do calculating

[traces,TT]=PassiveFW(model,source,SR,flag);

%% compare with theoretical green function in half space
source_id=1;
N_offset=10;
fm=TT.fm(source_id);
rr=TT.offset(source_id,:);

tmax=TT.time(end);
data1=TT.r(:,N_offset,source_id);
data2=TT.t(:,N_offset,source_id);
data3=TT.z(:,N_offset,source_id);
data=[data1 data2 data3];
[T, Urx, Utx, Uzz, Urz] = Lamb_3D(model.vs,model.vp,rr(N_offset),tmax);  %% step function response

sx=SR.sx(source_id);
sy=SR.sy(source_id);
rx=SR.rx(N_offset);
ry=SR.ry(N_offset);
[Azi,~]=cart2pol(rx-sx,ry-sy);
% sum displacements over ur,ut,uz in cylindrical coordinates

u1=(source.Fx(source_id)*cos(Azi)+source.Fy(source_id)*sin(Azi)).*Urx+source.Fz(source_id)*Urz;
u2=(source.Fy(source_id)*cos(Azi)-source.Fx(source_id)*sin(Azi)).*Utx;
u3=(source.Fx(source_id)*cos(Azi)+source.Fy(source_id)*sin(Azi)).*-Urz+source.Fz(source_id)*Uzz;

UU=[u1' u2' u3'];

for ii=1:3
[T1,zz_delta1]=Lamb_3D_compare(T,UU(:,ii),tmax,fm); %% Take the derivative
figure
plot(T1,zz_delta1/max(abs(zz_delta1)),'b-');
hold on;
plot(TT.time,data(:,ii)/max(abs(data(:,ii))),'r*');
legend_FontSize = legend('Analytical','Numerical');
end

toc

