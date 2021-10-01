function SR=S_R_geometry(source_number);
%**************************************************************************
%%%%%%%%%%%%  set ambient source-receiver geometry    %%%%%%%%%%%%%%%%%%%%%
%% traffic noise model
source_rang=300; %% radius of sources
rec_rang=100;    %% radius of spread
dx=5;  %% receiver spacing
sourceX= (2*(source_rang)*rand(1,source_number)-source_rang);  %% Source random  distribution in x-axis 
offline=10;
sourceY=1.5*offline.*rand(1,source_number)+offline/2;  %% offline= Line distance between road center line and receiving line
% sourceY=offline*ones(1,source_number);
receiverX=-rec_rang:dx:rec_rang;  %% receiver spread
receiverY=zeros(1,length(receiverX));  %% 1D spread

%**************************************************************************
%% love model
% source_rang=300; %% radius of sources
% rec_rang=100;    %% radius of spread
% dx=5;  %% receiver spacing
% sourceX= -150*ones(1,source_number); %% Source random  distribution in x-axis 
% sourceY=(2*(source_rang)*rand(1,source_number)-source_rang);
% receiverX=-rec_rang:dx:rec_rang;  %% receiver spread
% receiverY=zeros(1,length(receiverX));  %% 1D spread

%**************************************************************************
%%% model of uniform distribution
% source_rang1=500; %% radius of sources
% source_rang2=1000; %% radius of sources
% rec_rang=100;    %% radius of spread
% dx=5;  %% receiver spacing
% r=(source_rang2-source_rang1)*rand(1,source_number)+source_rang1;
% thR=2*pi*rand(1,source_number);
% [sourceX,sourceY]=pol2cart(thR,r);
% receiverX=-rec_rang:dx:rec_rang;  %% receiver spread
% receiverY=zeros(1,length(receiverX));  %% 1D spread
%**************************************************************************
SR.sx=sourceX;
SR.sy=sourceY;
SR.rx=receiverX;
SR.ry=receiverY;