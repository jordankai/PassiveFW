clear;
tic
%%  half-space model
% Damp_p=1./(2*Qp); % damping
% Damp_s=1./(2*Qs);

% model.vs   = 440;   % Shear wave velocity [m/s]
% model.vp   =760;
% model.dns  = 2;    % Mass density of layers (kg/m^3)
% model.thk  = inf;      % depth of layers [m]
% model.Damp_s=0;  %  viscoealstic damping
% model.Damp_p=0;  % A zero represents elastic media.

%% two-layer model
model.vs   = [220 440];   % Shear wave velocity [m/s]
model.vp   =[380 760];
model.dns  = [1 1 ]*2000;    % Mass density of layers (kg/m^3)
model.thk    = [10 inf];      % depth of layers [m]
Qs=[30 30];
Qp=[30 30];
model.Damp_p =[0 0] ;%1./(2*Qp); % damping
model.Damp_s = [0 0];%1./(2*Qs);

%%  three-layer model
% model.vs   = [300 180 450];   % Shear wave velocity [m/s]
% model.vp   =[600 1700 2000];
% model.dns  = [2 2 2 ];    % Mass density of layers (kg/m^3)
% model.thk  = [20 20 inf];      % depth of layers [m]
% model.Damp_s=[0 0 0];  %  viscoealstic damping
% model.Damp_p=[0 0 0];  % A zero represents elastic media.

%%  five-layer model   %  mi binbin
% vs = [270 367 125 453 540];   % Shear wave velocity [m/s]
% vp =[750 1400 550 1600 1800];
% dns = [1.86 1.91 1.96 2.02 2.09];    % Mass density of layers (kg/m^3)
% thk = [3 3 4 4 inf];      % depth of layers [m]
% Damp_s =[0 0 0 0 0];  %  viscoealstic damping
% Damp_p =[0 0 0 0 0];  % A zero represents elastic media.

source.fm=20;  %% source main frequency
source.dt=0.001;  %% time step
source.offset=1:200;  %% receivers position
source.maxfre=100;  %% max frequency to be calculated

flag='hv';
[ut,tt]=active_modeling_DWM(model,source,flag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot shot gathers
flag2=5;
wigb_normalized(ut(:,:,flag2),2,source.offset,tt);

%% dispersion imaging

 vv=100:500;
 [p1, f1]=phase_kai(ut(:,:,5),vv,source.offset,source.dt,50,5,'hankel');
 hold on
 plot(vr(:,1),vr(:,2:end),'*');
 
 vv=100:500;
 [p2, f2]=phase_kai(ut(:,:,2),vv,source.offset,source.dt,50,5,'hankel');
 hold on
 plot(vl(:,1),vl(:,2:end),'*');

toc
