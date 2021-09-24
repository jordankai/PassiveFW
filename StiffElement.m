function [Ksvp, Ksh] = StiffElement(Cp, Cs, Dp, Ds, rho, h, kr, w)
% Computes the SH and SVP stiffness matrices for a homogeneous layer
% Same sign convention as in Kausel, Elastodynamics
k = abs(kr);
wr = real(w);
wi = imag(w);
if wi>0, wi=-wi; end
w = abs(wr)+1i*wi;
% Cp = Cp*sqrt(1+2*1i*Dp*sign(wr));  
% Cs = Cs*sqrt(1+2*1i*Ds*sign(wr));

Cp = Cp*(1+1i*Dp*sign(wr));  
Cs = Cs*(1+1i*Ds*sign(wr));

a2 = (Cs/Cp)^2;
G = Cs^2*rho;
kp = sqrt(k^2-(w/Cp)^2);
if real(kp)<0
  kp = -kp;
end
ks = sqrt(k^2-(w/Cs)^2);
if real(ks)<0
  ks = -ks;
end
if k~=0
  p = sqrt(1-(w /(k*Cp))^2);
  s = sqrt(1-(w /(k*Cs))^2);
  if real(p)<0
    p = -p;
  end
  if real(s)<0
    s = -s;
  end
else
  p = 0;
  s = 0;
end
if h == inf
  if k>0 && w~=0
    Ksh = ks*G;
    Ksvp = 2*k*G*((1-s^2)/(2*(1-p*s))*[p -1;-1 s] + [0 1; 1 0]);
  elseif k>0 && w==0
    Ksh = k*G;
    Ksvp = 2*k*G/(1+a2)*[1  a2; a2  1];
  elseif k==0 && w~=0
    Ksh = 1i*w*rho*Cs;
    Ksvp = 1i*w*rho*[Cs  0; 0  Cp];
  elseif k==0 && w==0
    Ksh = 0;
    Ksvp = [0 0; 0 0];
  end
else
  if k>0 && w~=0
    RealPart_s = real(ks*h);
    ImagPart_s = imag(ks*h);
    RealPart_p = real(kp*h);
    ImagPart_p = imag(kp*h);
    if RealPart_s >= 0
      Shs = ((1-exp(-2*RealPart_s))*cos(ImagPart_s) + (1+exp(-2*RealPart_s))*sin(ImagPart_s)*1i)/2;
      Chs = ((1+exp(-2*RealPart_s))*cos(ImagPart_s) + (1-exp(-2*RealPart_s))*sin(ImagPart_s)*1i)/2;
    else
      Shs = -((1-exp(2*RealPart_s))*cos(ImagPart_s) - (1+exp(2*RealPart_s))*sin(ImagPart_s)*1i)/2;
      Chs =  ((1+exp(2*RealPart_s))*cos(ImagPart_s) - (1-exp(2*RealPart_s))*sin(ImagPart_s)*1i)/2;
    end
    if RealPart_p >= 0
      Shp = ((1-exp(-2*RealPart_p))*cos(ImagPart_p) + (1+exp(-2*RealPart_p))*sin(ImagPart_p)*1i)/2;
      Chp = ((1+exp(-2*RealPart_p))*cos(ImagPart_p) + (1-exp(-2*RealPart_p))*sin(ImagPart_p)*1i)/2;
    else
      Shp = -((1-exp(2*RealPart_p))*cos(ImagPart_p) - (1+exp(2*RealPart_p))*sin(ImagPart_p)*1i)/2;
      Chp =  ((1+exp(2*RealPart_p))*cos(ImagPart_p) - (1-exp(2*RealPart_p))*sin(ImagPart_p)*1i)/2;
    end
    Cths = Chs / Shs;
    ABS_Real_s = abs(real(RealPart_s));
    ABS_Real_p = abs(real(RealPart_p));
    kh = k*h;
    ps = p*s;
    if p == 0
      Denom = 2*(exp(-ABS_Real_s)-Chs)+Shs*k*h/s;
    elseif s == 0
      Denom = 2*(exp(-ABS_Real_p)-Chp)+Shp*k*h/p;
    else
      Denom = 2*(exp(-(ABS_Real_s+ABS_Real_p))-Chp*Chs)+(1/ps+ps)*Shp*Shs;
    end
    Chp_Chs_D = (Chp*Chs)/Denom;
    Chp_Shs_D = (Chp*Shs)/Denom;
    Shp_Chs_D = (Shp*Chs)/Denom;
    Shp_Shs_D = (Shp*Shs)/Denom;
    Chp_D = exp(-ABS_Real_s)*Chp/Denom;
    Shp_D = exp(-ABS_Real_s)*Shp/Denom;
    Chs_D = exp(-ABS_Real_p)*Chs/Denom;
    Shs_D = exp(-ABS_Real_p)*Shs/Denom;
    Um_D = exp(-ABS_Real_p-ABS_Real_s)/Denom;
    Shs = Shs*exp(ABS_Real_s);
    Ksh = ks*G*[Cths -1/Shs; -1/Shs Cths];
    S1 = k*G*(1-s^2);
    S2 = k*G*(1+s^2);
    if p == 0
      k11 = S1*(Chp_Shs_D/s);
      k12 = S1*(Um_D-Chp_Chs_D) + S2;
      k22 = S1*(kh*Chs_D -s*Chp_Shs_D);
    elseif s == 0
      k11 = S1*(kh*Chp_D-p*Shp_Chs_D);
      k12 = S1*(Um_D-Chp_Chs_D) + S2;
      k22 = S1*(Shp_Chs_D/p);
    else
      k11 = S1*(Chp_Shs_D - ps*Shp_Chs_D)/s;
      k12 = S1*(Um_D-Chp_Chs_D+ps*Shp_Shs_D)+ S2;
      k22 = S1*(Shp_Chs_D - ps*Chp_Shs_D)/p; 
    end
    K11 = [k11,  k12;  k12, k22];
    K22 = [k11, -k12; -k12, k22];
    if p == 0
      K12 = S1* [-Shs_D/s,          Chp_D-Chs_D; -(Chp_D-Chs_D),  s*Shs_D-kh*Um_D];
    elseif s == 0
      K12 = S1*[ p*Shp_D-kh*Um_D,   Chp_D-Chs_D; -(Chp_D-Chs_D),  -Shp_D/p];
    else
      K12 = S1*[(ps*Shp_D-Shs_D)/s, Chp_D-Chs_D; -(Chp_D-Chs_D),  (ps*Shs_D-Shp_D)/p];
    end
    K21 = [K12(1,1) K12(2,1); K12(1,2) K12(2,2)];
    Ksvp = [K11 K12; K21 K22];
  elseif k>0 && w==0
    kh = k*h;
    RealPart = real(kh);
    ImagPart = imag(kh);
    if RealPart >= 0
      Sh = ((1-exp(-2*RealPart))*cos(ImagPart) + (1+exp(-2*RealPart))*sin(ImagPart)*1i)/2;
      Ch = ((1+exp(-2*RealPart))*cos(ImagPart) + (1-exp(-2*RealPart))*sin(ImagPart)*1i)/2;
    else
      Sh = -((1-exp(2*RealPart))*cos(ImagPart) - (1+exp(2*RealPart))*sin(ImagPart)*1i)/2;
      Ch =  ((1+exp(2*RealPart))*cos(ImagPart) - (1-exp(2*RealPart))*sin(ImagPart)*1i)/2;
    end
    Cth = Ch / Sh;
    Denom = (1+a2)^2*Sh^2 - kh^2*(1-a2)^2*exp(-2*abs(RealPart));
    Ch_Ch_D = Ch*Ch / Denom;
    Sh_Sh_D = Sh*Sh / Denom;
    Ch_Sh_D = Ch*Sh / Denom;
    Ch_D = Ch*exp(-abs(RealPart)) /Denom;
    Sh_D = Sh*exp(-abs(RealPart)) /Denom;
    Um_D =  exp(-2*abs(RealPart)) /Denom;
    Ch = Ch*exp(abs(RealPart));
    Sh = Sh*exp(abs(RealPart));
    Ksh = k*G*[Cth -1/Sh; -1/Sh Cth];
    K11 = 2*k*G*([(1+a2)*Ch_Sh_D-kh*(1-a2)*Um_D,  -(1+a2)*Sh_Sh_D; ...
      -(1+a2)*Sh_Sh_D,   (1+a2)*Ch_Sh_D+kh*(1-a2)*Um_D] + [0 1; 1 0]);
    K22 = [K11(1,1), -K11(1,2); -K11(2,1), K11(2,2)];
    K12 = 2*k*G*[kh*(1-a2)*Ch_D-(1+a2)*Sh_D kh*(1-a2)*Sh_D;  -kh*(1-a2)*Sh_D -(kh*(1-a2)*Ch_D+(1+a2)*Sh_D)];
    K21 = [K12(1,1) K12(2,1); K12(1,2) K12(2,2)];
    Ksvp = [K11 K12; K21 K22];
  elseif k==0 && w~=0    % elseif k==0 && w~=0
    omegaP = w*h/Cp; omegaS = w*h/Cs;
    Ksh = rho*Cs*w*[cot(omegaS) -1/sin(omegaS); -1/sin(omegaS) cot(omegaS)];
    Ksvp = rho*w*[Cs*cot(omegaS), 0, -Cs/sin(omegaS), 0; ...
                   0, Cp*cot(omegaP), 0, -Cp/sin(omegaP); ...
                   -Cs/sin(omegaS), 0, Cs*cot(omegaS), 0; ...
                   0, -Cp/sin(omegaP), 0, Cp*cot(omegaP)];
  elseif k==0 && w==0    % elseif k==0 && w==0
    a_2 = 1/a2;
    Ksh = G/h*[1 -1; -1 1];
    Ksvp = G/h*[1 0 -1 0; 0 a_2 0 -a_2; -1 0 1 0; 0 -a_2 0 a_2];
  end
end
% Modify matrices for negative wavenumber and/or frequency
if kr<0
  n = length(Ksvp);
  if n>2
    Ksh(1,2)  = -Ksh(1,2);
    Ksh(2,1)  = -Ksh(2,1);
  end
  for m=2:2:n
    Ksvp(m,:) = -Ksvp(m,:);
  end
  for m=2:2:n
    Ksvp(:,m) = -Ksvp(:,m);
  end
end
if wr<0
  Ksh  = conj(Ksh);
  Ksvp = conj(Ksvp);
end
end