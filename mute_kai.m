function ut=mute_kai(ut,c,x,dt);

index_kk=fix(x/c/dt);
for ii=1:length(x)
    ut(1:index_kk(ii),ii)=0;
end