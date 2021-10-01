
function  ut2=conv_wavelet(t,ut,fm)
% t: time vector row
% fm: main frequency row
% ut: signal, one colum represents one trace
[m,n]=size(ut);
if length(fm)==1
   fm=ones(1,n)*fm; 
end
dt=t(2)-t(1);
t5=t(1):dt:(m-1)*dt;
ss3=zeros(length(t5),n);
ss3(1:length(t),:)=ut;
pz2=rickerkai(t5,fm);
ut=real(ifft(fft(ss3).*fft(pz2)));
ut2=ut(1:length(t),:);




