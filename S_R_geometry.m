function SR=S_R_geometry(source_number);
%**************************************************************************
%%%%%%%%%%%%  set ambient source-receiver geometry    %%%%%%%%%%%%%%%%%%%%%
source_rang=300; %% radius of sources
rec_rang=100;    %% radius of spread
dx=5;  %% receiver spacing
sourceX= -150*ones(1,source_number);%(2*(source_rang)*rand(1,source_number)-source_rang);  %% Source random  distribution in x-axis %
% sourceX(1,1:ceil(source_number/2))=(source_rang-rec_rang)+rec_rang*rand(1,ceil(source_number/2));
% sourceX(1,ceil(source_number/2)+1:source_number)=-((source_rang-rec_rang)+rec_rang*rand(1,ceil(source_number/2)));
offline=10;
sourceY=(2*(source_rang)*rand(1,source_number)-source_rang);%2*offline.*rand(1,source_number);  %% offline= Line distance between road center line and receiving line
%sourceY=offline*ones(1,source_number);
receiverX=-rec_rang:dx:rec_rang;  %% receiver spread
receiverY=zeros(1,length(receiverX));  %% 1D spread

%**************************************************************************
%%% model pang ,uniform distribution
% source_rang1=500; %% radius of sources
% source_rang2=1000; %% radius of sources
% rec_rang=100;    %% radius of spread
% dx=5;  %% receiver spacing
% r=(source_rang2-source_rang1)*rand(1,source_number)+source_rang1;
% thR=2*pi*rand(1,source_number);
% [sourceX,sourceY]=pol2cart(thR,r);
% % offline= Line distance between road center line and receiving line
% receiverX=-rec_rang:dx:rec_rang;  %% receiver spread
% receiverY=zeros(1,length(receiverX));  %% 1D spread
%**************************************************************************
SR.sx=sourceX;
SR.sy=sourceY;
SR.rx=receiverX;
SR.ry=receiverY;