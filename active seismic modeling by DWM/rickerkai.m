function pz = rickerkai(t,fm)
% t: time vector row
% fm: main frequency row

%% max value at t0
t0=2./fm;
pz=(1-2*(pi*fm.*(t'-t0)).^2).*exp(-(pi*fm.*(t'-t0)).^2);

end