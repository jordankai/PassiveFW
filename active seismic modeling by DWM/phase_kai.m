function [p1, f1]=phase_kai(data_kai,v,x,dt,f_max,varargin)
%% Dispersion imaging by the Frequency‐Bessel or Hankel transform or phase-shift
%**************************************************************************
% Input arguments
%  data_kai   = seismic traces: each column represents one traces
%  v = a vector of velocity, one row
%  x = a vector of offset, one row
%  dt = time step
%  f_max = max. frequency to be calculated
%  df_factor  = Magnification time duration, so df is more finer
%  flag  = 'bessel', 'hankel' or 'phase'

%**************************************************************************
% output arguments
% p1 = Dispersion image matrix with velocity × frequency
% f1 = Frequency row

% example: [p1 f1]=phase_kai(ut5,100:600,1:200,0.001,100,1,'hankel','z');
%**************************************************************************

if  nargin==5
    df_factor=1;
    flag='phase';
    flag_component='V';
end
if  nargin==6
    df_factor=varargin{1};
    flag='phase';
    flag_component='V';
end
if nargin==7
    df_factor=varargin{1};
    flag=varargin{2};
    flag_component='V';
elseif nargin>7
    df_factor=varargin{1};
    flag=varargin{2};
    flag_component=varargin{3};
end

if ~any(strcmp(flag,{'hankel','bessel','phase'}))
    error('Wrong input strings within ({hankel,bessel,phase})');
end
if ~any(strcmp(flag_component,{'v','V','r','R','t','T'}))
    error('Wrong input strings within ({v,r,t,V,R,T})');
end

time_num=length(data_kai(:,1))*df_factor;
df=1/(dt*time_num);
f=df:df:ceil(f_max/df)*df;
f1=0:df:ceil(f_max/df)*df;
n=length(data_kai(1,:));
fn=length(f);
vn=length(v);
p1=zeros(vn,fn+1);

seis_factor=zeros(time_num,n);
seis_factor(1:length(data_kai(:,1)),:)=data_kai;
for i=1:n
    if strcmp(flag,{'bessel'}) || strcmp(flag,{'hankel'})
        a(i,:)=(fft(seis_factor(:,i)));
    else
        a(i,:)=angle(fft(seis_factor(:,i)));
    end
end

w=2*pi*f;
p=zeros(vn,fn);
%%**************************************************************************
% %Trapezoidal integral
% for k=1:n
%     be1(:,:,k)=besselj(0,w*x(k)./v')*x(k).*a(k,2:(fn+1)); % Bessel transform method
%     be1(:,:,k)=besselh(0,1,w*x(k)./v')*x(k).*a(k,2:(fn+1)); % Hankel transform method
% end
% p=trapz(x(1:n),be1,3);
%%**************************************************************************
%% Ref: The Frequency‐Bessel Spectrograms of Multicomponent
%% Cross-Correlation Functions From Seismic
%% Ambient Noise
%% fast method
rrspan=[0 x x(end)];
kk=w./v';
for ii=2:length(x)+1
    % ref: Inversion of shallow-seismic wavefields: I.Wavefield transformation
    % Modified frequency–Bessel transform method for dispersion imaging of Rayleigh waves from ambient seismic noise
    if strcmp(flag,{'bessel'}) % Bessel transform method
        newr2=(rrspan(ii+1)^2+2*rrspan(ii)*(rrspan(ii+1)-rrspan(ii-1))-rrspan(ii-1)^2)/8;
        hh=(a(ii-1,2:(fn+1))); %% because matlab fft is exp(-iwt),while exp(iwt) in ref
        if  any(strcmp(flag_component,{'r','R'}))
            IIw=newr2*hh.*besselj(1,kk*rrspan(ii));
        else
            IIw=newr2*hh.*besselj(0,kk*rrspan(ii));
        end
        p=p+IIw;
    elseif strcmp(flag,{'hankel'})  % hankel transform method
        newr2=(rrspan(ii+1)^2+2*rrspan(ii)*(rrspan(ii+1)-rrspan(ii-1))-rrspan(ii-1)^2)/8;
        hh=(a(ii-1,2:(fn+1)));
        if any(strcmp(flag_component,{'r','R'}))
            IIw=newr2*(hh.*besselh(1,1,kk*rrspan(ii)));
            %IIw=newr2*hh.*0.5.*(w.^2).*besselh(1,2,kk*rrspan(ii));
        else
            IIw=newr2*(hh.*besselh(0,1,kk*rrspan(ii)));
        end
        p=p+IIw;
    else
        p=p+exp(sqrt(-1)*(w*x(ii-1)./v'+a(ii-1,2:fn+1)));  % phase-shift method
    end
    
end
% normalization in frequency
for ii=1:fn
    if max(abs(p(:,ii)))~=0
        p1(:,ii+1)=abs(p(:,ii))/max(abs(p(:,ii)));
    end
end
% plot
figure;
colormap(jet);
imagesc(f1,v,p1)
set(gca,'YDir','normal')




