

function [VV_spectra,VR_spectra,HR_spectra,HT_spectra,HV_spectra]=Green_spectra_k(Cs,Cp,Rho,H,Damp_s,Damp_p,kk,ff);

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
nL = length(Cs);             % Number of layers (materials)
nL2 = 2*nL;                  % double the number of layers
qR = zeros(nL2,1);  qR(1)=1;        % load vector for Radial
qT = zeros(nL,1);   qT(1)=1;        % load vector for Tangential
qV = zeros(nL2,1);  qV(2)=1;        % load vector for Vertical
VV_spectra=zeros(length(kk),length(ff));
VR_spectra=zeros(length(kk),length(ff));
HV_spectra=zeros(length(kk),length(ff));
HR_spectra=zeros(length(kk),length(ff));
HT_spectra=zeros(length(kk),length(ff));
parfor ii=1:length(ff)
    w_ew=2*pi*ff(ii);    
    %**************************************************************************
    [KGsh, KGsvp]=arrayfun(@(ky)(GlobalStiffMat(Cp,Cs,Damp_p,Damp_s,Rho,H,ky,w_ew)),kk,'UniformOutput',false);
    UA_temp_z=cellfun(@(x) x\qV,KGsvp,'UniformOutput',false);
    %   UA_temp=cellfun(@(x) pinv(x)*qR,KGsvp,'UniformOutput',false);%% pinv is not recommended,time consuming
    UA_z=cell2mat(UA_temp_z);
    
    if ~any(isnan(UA_z))
        VR_spectra(:,ii) = (UA_z(1,:));
        VV_spectra(:,ii) = (UA_z(2,:));
    end
    
    UA_temp_z=cellfun(@(x) x\qR,KGsvp,'UniformOutput',false);
    
    UA_z=cell2mat(UA_temp_z);
    
    UA_temp_t=cellfun(@(x) x\qT,KGsh,'UniformOutput',false);
    UA_t=cell2mat(UA_temp_t);
    
    if ~any(isnan(UA_z))
        HR_spectra(:,ii) = (UA_z(1,:));
        HV_spectra(:,ii) = (UA_z(2,:));
    end
    if ~any(isnan(UA_t))
        HT_spectra(:,ii) = (UA_t(1,:));
    end
    
    
    
end

% VR_spectra=VR_spectra./max(VR_spectra);
% VV_spectra=VV_spectra./max(VV_spectra);
% HR_spectra=HR_spectra./max(HR_spectra);
% HT_spectra=HT_spectra./max(HT_spectra);
% HV_spectra=HV_spectra./max(HV_spectra);








