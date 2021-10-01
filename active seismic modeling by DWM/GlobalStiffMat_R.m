function KGsvp = GlobalStiffMat_R(Cp, Cs, Dp, Ds, rho, h, k, w)
% Assembles the stiffness matrix for the layered soil
%% Check properties
Cp = Cp(:);
Cs = Cs(:);
Dp = Dp(:);
Ds = Ds(:);
rho = rho(:);
h = h(:);
error = 0;
nCp = length(Cp);
nCs = length(Cs);
if nCp ~= nCs
    error = 1;
end
nDp = length(Dp);
if nCp ~= nDp
    error = 1;
end
nDs = length(Ds);
if nCp ~= nDs
    error = 1;
end
nrho = length(rho);
if nCp ~= nrho
    error = 1;
end
nh = length(h);
if nCp ~= nh
    error = 1;
end
if error
    disp('Error: properties dimension mismatch');

    KGsvp = [];
    return;
end

%% Formation of matrices

KGsvp = zeros(2*nCp, 2*nCp);

%% Assemblage of top layers
for iLayer = 1:nCp-1
    Ksvp = StiffElement_R(Cp(iLayer), Cs(iLayer), Dp(iLayer), Ds(iLayer), rho(iLayer), h(iLayer), k, w);
  
    KGsvp((iLayer-1)*2+1:(iLayer+1)*2, (iLayer-1)*2+1:(iLayer+1)*2) = KGsvp((iLayer-1)*2+1:(iLayer+1)*2, (iLayer-1)*2+1:(iLayer+1)*2) + Ksvp;
end

%% Assemblage of last element
Ksvp = StiffElement_R(Cp(nCp), Cs(nCp), Dp(nCp), Ds(nCp), rho(nCp), h(nCp), k, w);

KGsvp(2*nCp-1:2*nCp, 2*nCp-1:2*nCp) = KGsvp(2*nCp-1:2*nCp, 2*nCp-1:2*nCp) + Ksvp(1:2, 1:2);
end

