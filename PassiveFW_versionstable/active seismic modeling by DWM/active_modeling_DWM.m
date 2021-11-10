
function  [ut,tt]=active_modeling_DWM(model,source,flag);

% The code is based on E.Kausel's moving bell load code
%**************************************************************************
%    (many thanks to Professor E.Kausel)
%    Written and Copyleft by Eduardo Kausel, MIT, kausel@mit.edu
%    Version 1, May 21, 2017
%    You may use and modify this code, provided you acknowledge the source
%    Referrence "Generalized stiffness matrix method for layered soils" (2018)
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
%  model   = Material properties
%  source  = source and geometry properties
%  flag  = 'V' or 'v' for vertical load and 'H' or 'h' for horizontal load
%          'all' or 'All' for both of them
%**************************************************************************
% output arguments
% ut(:,:,1)-Radial compenent due to Horizontal point load
% ut(:,:,2)-Tangential compenent due to Horizontal point load
% ut(:,:,3)-Vertical compenent due to Horizontal point load
% ut(:,:,4)-Tangential compenent due to vertical point load
% ut(:,:,5)-Vertical compenent due to vertical point load
% tt-time sequence
%**************************************************************************

fm=source.fm;  %% source main frequency
dt=source.dt;  %% time step
offset=source.offset; %% receivers position
maxfre=source.maxfre; %% max frequency to be calculated

t0=2/fm;  %% source time delay 
delay_t=2*t0;  %% time delay of green function 
NN_green=fix(delay_t/dt);  

Td=2*max(offset)/min(model.vs)/0.862+t0+delay_t; %% total record time
Nt=ceil(Td/dt); 

%% calculation of green function in f-k domain
UK=Green_HV_delta_one_par(model,Nt,dt,max(offset),delay_t,maxfre,flag);

%% do the Hankel transform from f-k domain to f-space domain using discrete wavenumber method 
[Wavef,tt2]=Green_HV_delta_two_par(UK,Nt,dt,offset,maxfre,flag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convolve with Ricker wavelet
tt=tt2(1:end-NN_green);
ut=zeros(length(tt),length(offset),5);
for ii=1:5
ut_temp=conv_wavelet(tt2,Wavef(:,:,ii),fm);
ut(:,:,ii)=ut_temp(NN_green+1:end,:); %% remove time delay of Green function 
end





