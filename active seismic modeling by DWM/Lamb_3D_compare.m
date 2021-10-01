%% compare with analytic solution (Kausel 2012)

function  [T1,zz_delta1]=Lamb_3D_compare(Tr,Urz,tmax,fm)

index=find(abs(Urz)==inf);  %% delete inf value for step response
Urz(index)=[];
Tr(index)=[];

index1=isnan(Urz);
Urz(index1)=[];
Tr(index1)=[];

%% method 1
T1=0:min(diff(Tr)):(tmax-min(diff(Tr)));
%T1=0:(t(2)-t(1)):(tmax-min(diff(Tr)));
ss1=interp1(Tr,Urz,T1);
for ii=2:length(ss1)-1
    ss2(ii)=(ss1(ii+1)-ss1(ii-1))/2;
end
ss2(1)=0;
ss2(length(ss1))=0;

pz1=rickerkai(T1,fm);
zz_delta1=real(ifft(fft(ss2).*fft(pz1')));




