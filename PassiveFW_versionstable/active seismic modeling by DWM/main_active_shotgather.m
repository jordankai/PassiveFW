clear;
tic
%%  half-space model
% Damp_p=1./(2*Qp); % damping
% Damp_s=1./(2*Qs);

model.vs   = 440;   % Shear wave velocity [m/s]
model.vp   =760;
model.dns  = 2;    % Mass density of layers (kg/m^3)
model.thk  = inf;      % depth of layers [m]
model.Damp_s=0;  %  viscoealstic damping
model.Damp_p=0;  % A zero represents elastic media.

%% two-layer model
% model.vs   = [220 440];   % Shear wave velocity [m/s]
% model.vp   =[380 760];
% model.dns  = [1 1 ]*2000;    % Mass density of layers (kg/m^3)
% model.thk    = [10 inf];      % depth of layers [m]
% Qs=[30 30];
% Qp=[30 30];
% model.Damp_p =[0 0] ;%1./(2*Qp); % damping
% model.Damp_s = [0 0];%1./(2*Qs);

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

%% compare with theoretical green function in half space
N_offset=10;
r=source.offset;
[T, Urx, Utx, Uzz, Urz] = Lamb_3D(model.vs,model.vp,r(N_offset),tt(end));  %% step function response
UU=[Urx' Utx' -Urz' Urz' Uzz'];
for ii=1:5
[T1,zz_delta1]=Lamb_3D_compare(T,UU(:,ii),tt(end),source.fm); %% Take the derivative
figure
plot(T1,zz_delta1/max(abs(zz_delta1)),'b-');
hold on;
data=ut(:,N_offset,ii);
plot(tt,data/max(abs(data)),'r*');
legend_FontSize = legend('Analytical','Numerical');
end


% FigFontSize=10.5;
% FigWidth=3.33;   FigHeight=2.3;
% figure
% plot(T1,zz_delta1/max(abs(zz_delta1)),'b-','Linewidth',1);
% hold on;
% plot(tt,data/max(abs(data)),'r*','MarkerIndices',[1:15:91 93:3:165 185:20:500],'MarkerSize',4);
% 
% legend_FontSize = legend('Analytical','Numerical');
% set(legend_FontSize,'FontSize',9) 
% set(gca, 'FontSize', FigFontSize)
% set(gca,'LineWidth',1)
% xlabel('Time (s)','Fontname', ' Arial ','FontSize',FigFontSize)
% set(gca,'XLim',[0 0.5])
% % set(gca,'XTick',[0:0.1:0.5])
% ylabel('Normalized Amplitude','Fontname', ' Arial ','FontSize',FigFontSize)
% % set(gca,'YTick',[0.2:0.2:1])
% set(gcf,'units','inch')
% pos = [5, 2, FigWidth, FigHeight];
% set(gcf,'position',pos);
% set(gca,'Position',[.17 .2 0.78 0.78]);

toc
