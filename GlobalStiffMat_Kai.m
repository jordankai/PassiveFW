
function UU=GlobalStiffMat_Kai(Cp,Cs,Damp_p,Damp_s,Rho,H,ky,w_ew)

[KGsh, KGsvp]=arrayfun(@(w)(GlobalStiffMat(Cp,Cs,Damp_p,Damp_s,Rho,H,ky,w)),w_ew,'UniformOutput',false);
qR=zeros(2*length(Cs),1);  qR(1,1)=1;
qT=zeros(length(Cs),1);    qT(1,1)=1;
for ii=1:length(w_ew)
    UA_temp=KGsvp{ii}\qR;
    UR(ii) = UA_temp(1);
    
    UA_temp=KGsh{ii}\qT;
    UT(ii) = UA_temp(1);
end
UU=UR-UT;
