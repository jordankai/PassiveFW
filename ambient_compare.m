
flag=3;
Cp=760;
Cs=440;
sources_number=2;
rr=offset(sources_number,:);
fm=10;

N_offset=11;
tmax=t_one(end);
data=TT_z(:,N_offset,sources_number);
[T, Urx, Utx, Uzz, Urz] = Lamb_3D(Cs,Cp,rr(N_offset),tmax);  %% step function response
%Lamb_3D_compare(data,t_one,T,Urz,tmax,fm,flag);

[tha,~]=cart2pol(-50-100,0-100)
u1=(1*cos(tha)+0*sin(tha))*Urx+0*Urz;
u2=(-1*sin(tha)+0*cos(tha))*Utx;
u3=(1*cos(tha)+0*sin(tha))*(-1)*Urz+0*Uzz;
Lamb_3D_compare(data,t_one,T,u3,tmax,fm,flag);