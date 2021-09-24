
%% model setting
model.Cs   = [220 440];   % Shear wave velocity [m/s]
model.Cp   =[380 760];
model.Rho  = [1 1 ]*2000;    % Mass density of layers (kg/m^3)
model.H    = [10 inf];      % depth of layers [m]
model.Qs=[30 30];
model.Qp=[30 30];
model.Damp_p =[0 0] ;%1./(2*Qp); % damping
model.Damp_s = [0 0];%1./(2*Qs);

%% time and source setting
source_number=1000;
delka_t=1000; % record time(s)
dt=0.005;  
fmax=50;
fm1=5;  %% min. Dominant Frequency of Ricker wavelet
fm2=30;  %% max. Dominant Frequency of Ricker wavelet

%% Geometry setting
SR=S_R_geometry(source_number);

%% do calculating
[traces,TT]=PassiveFW_all(dt,fmax,fm1,fm2,delka_t,source_number,SR,model);

