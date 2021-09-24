%clear;
%%%%  three-layer model
vs3   = [300 180 450];   % Shear wave velocity [m/s]
vp3   =[600 1700 2000];
dns3  = [2 2 2 ];    % Mass density of layers (kg/m^3)
thk3    = [20 20 inf];      % depth of layers [m]
Damp_s3=[0 0 0];  %  viscoealstic damping
Damp_p3=[0 0 0];  % A zero represents elastic media.

%%%%%%%%  five-layer model   %  mi binbin
vs5   = [270 367 125 453 540];   % Shear wave velocity [m/s]
vp5   =[750 1400 550 1600 1800];
dns5  = [1.86 1.91 1.96 2.02 2.09];    % Mass density of layers (kg/m^3)
thk5    = [3 3 4 4 inf];      % depth of layers [m]
Damp_s5=[0 0 0 0 0];  %  viscoealstic damping
Damp_p5=[0 0 0 0 0];  % A zero represents elastic media.


fm=20;
t0=2/fm;
delay_t=2*t0;
dt=0.001;
NN=fix(delay_t/dt);

offset=1:200;
Td=2*max(offset)/min(vs5)/0.862+t0+delay_t;
Nt=ceil(Td/dt);

[Ky,uk1,uk2,uk3]=Green_HV_delta_one(vs5,vp5,dns5,thk5,Damp_s5,Damp_p5,Nt,dt,max(offset),delay_t,100,'h');
[Wavef1,tt]=Green_HV_delta_two(Nt,dt,Ky,uk1,uk2,uk3,offset,100,'h');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot v_v
flag=2;
ut5=conv_wavelet(tt,Wavef1(:,:,flag),fm);
ut5=ut5(NN+1:end,:);
tt5=tt(1:end-NN);

wigb_normalized(ut5,2,offset,tt5);


%% pick the dispersion peaks
vv=100:600;
[p1, f1]=phase_kai(ut5(2:2001,1:100),vv,1:100,0.001,70,5,'hankel');
nn=81;
dis_curve5(:,1)=f1(nn:end);
for ii=1:621
    [~,ind]=max(p1(:,nn+ii-1));
    dis_curve5(ii,2)=vv(1)+(ind-1)*(vv(2)-vv(1));
end
hold on;
plot(dis_curve5(:,1),dis_curve5(:,2),'*');

rms3=min(abs(dis_curve3(:,2)-cc_true3(:,2:8)),[],2);

