%%  half-space model
tic
% Damp_p=1./(2*Qp); % damping
% Damp_s=1./(2*Qs);

% model.vs   = 440;   % Shear wave velocity [m/s]
% model.vp   =760;
% model.dns  = 2;    % Mass density of layers (kg/m^3)
% model.thk  = inf;      % depth of layers [m]
% model.Damp_s=0;  %  viscoealstic damping
% model.Damp_p=0;  % A zero represents elastic media.

%%  three-layer model
% model.vs   = [300 180 450];   % Shear wave velocity [m/s]
% model.vp   =[600 1700 2000];
% model.dns  = [2 2 2 ];    % Mass density of layers (kg/m^3)
% model.thk  = [20 20 inf];      % depth of layers [m]
% model.Damp_s=[0 0 0];  %  viscoealstic damping
% model.Damp_p=[0 0 0];  % A zero represents elastic media.

%% two-layer model
model.vs   = [220 440];   % Shear wave velocity [m/s]
model.vp   =[380 760];
model.dns  = [1 1 ]*2000;    % Mass density of layers (kg/m^3)
model.thk    = [10 inf];      % depth of layers [m]
Qs=[30 30];
Qp=[30 30];
model.Damp_p =[0 0] ;%1./(2*Qp); % damping
model.Damp_s = [0 0];%1./(2*Qs);

%% time and source setting
source.number=1000; 
source.recordlength=1000; % record time(s)
source.dt=0.005; %  time step
source.maxfre=100;  %% max frequency to be calculated
source.fmin=5;  %% min. Dominant Frequency of Ricker wavelet
source.fmax=30;  %% max. Dominant Frequency of Ricker wavelet

%%%%%%%%%%%%  set point loads   %%%%%%%%%%%%%%%%%%%%%
flag='v';

if any(strcmp(flag,{'v','V'}))
    source.Fx=zeros(1,source.number); % point load in x-axis
    source.Fy=zeros(1,source.number); % point load in y-axis
    source.Fz=-rand(1,source.number);% point load in z-axis,z-axis is positive in the upward direction
elseif any(strcmp(flag,{'h','H'}))
    source.Fx=rand(1,source.number); % point load in x-axis
   % source.Fx=zeros(1,source.number); % point load in x-axis
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

toc

