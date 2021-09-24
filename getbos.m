function bos = getbos(vs,vp)
% get Poisson's ratio from Vp and Vs
temp = (vp./vs).^2;
bos = (1-0.5*temp)./(1-temp);
end