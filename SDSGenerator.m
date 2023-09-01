% Code for generating SDS Equations
% Yashica Khatri
% Setting up SDS eqn xalculator from Inkwan Park to get separated equations
% for short period variations

clc
clear 
% close all

syms l g h L G H ks Ks hs cs sSun Ls Gs MU
syms J2 b mu RE nu rss r rfun(L,G,H,l) crLpart crGpart crHpart crlpart TAfun(L,G, H,l) TA cfLpart cfGpart cflpart cfHpart Ecc Eccfun(L,G,H,l) cELpart cEGpart cElpart cEHpart E

getHtilde = 0;
getK = 0;
getKLW1 =  0;
getKLW2 = 0;
flagGetInitialDelOffset = 1;
flagGetFinalDelOffset = 1;
writeNewFile = 1;

x = [l;g;h;ks]; X = [L;G;H;Ks];
Xx = [X;x];
a = L^2/mu;
e = sqrt(1-(G/L)^2);
s = sqrt(1-(H/G)^2);
c = H/G;
n = mu^2/L^3;
eta = G/L;

drdM = a * e * sin(TA) / eta;
drda = r / a;
dadL = 2 / a / n;
drde = -a * cos(TA);
dedL = eta^2 / n / a^2 / e;
dMdL = 0; 
dedG = -eta / n / a^2 / e;
dMdG = 0; 
dMdH = 0; 
dEdM = a/r;
dEde = a/r * sin(Ecc);
dMdl = 1;
dfde = (2 + e * cos(TA)) * sin(TA) / eta^2;
dfdM = (a/r)^2 * eta;
dedl = 0;
dedH = 0;
dadl = 0;
dadG = 0;
dadH = 0;

crlpartVal = drde * dedl + drdM * dMdl; 
crLpartVal = drde * dedL + drdM * dMdL; 
crGpartVal = drde * dedG + drdM * dMdG; 
crHpartVal = drde * dedH + drdM * dMdH; 

cElpartVal = dEde * dedl + dEdM * dMdl;
cELpartVal = dEde * dedL + dEdM * dMdL;
cEGpartVal = dEde * dedG + dEdM * dMdG;
cEHpartVal = dEde * dedH + dEdM * dMdH;

cflpartVal = dfde * dedl + dfdM * dMdl;
cfLpartVal = dfde * dedL + dfdM * dMdL; 
cfGpartVal = dfde * dedG + dfdM * dMdG; 
cfHpartVal = dfde * dedH + dfdM * dMdH; 


Tc = 1/4*((c+1)*(cs+1)*cos(g+h-hs-ks)...
    +(c-1)*(cs-1)*cos(g-h+hs-ks)...
    -(c+1)*(cs-1)*cos(g+h-hs+ks)...
    -(c-1)*(cs+1)*cos(g-h+hs+ks)...
    +4*s*sSun*sin(g)*sin(ks));
Ts = 1/4*((-c-1)*(cs+1)*sin(g+h-hs-ks)...
    -(c-1)*(cs-1)*sin(g-h+hs-ks)...
    +(c+1)*(cs-1)*sin(g+h-hs+ks)...
    +(c-1)*(cs+1)*sin(g-h+hs+ks)...
    +4*s*sSun*cos(g)*sin(ks));

subPartialsInputs = [r, crLpart, crGpart, crHpart, crlpart, cfLpart, cfGpart, cfHpart, cflpart, TA, cELpart, cEGpart, cElpart, cEHpart, Ecc,crLpartVal,crGpartVal,crHpartVal,crlpartVal,cfLpartVal,cfGpartVal,cfHpartVal, cflpartVal,cELpartVal,cEGpartVal,cElpartVal,cEHpartVal];

%% OG Hamiltonian
H0 = -mu / 2 / a;
H2 = 2 / E^2 * nu * Ks;
H3const = -mu / 4 / a /eta^2 * RE^2 / r^2 * J2;
H3 = factorial(3) / E^3 * H3const * (2-3*s^2 + 3 * s^2 * cos(2*(TA+g)));
H5 = - factorial(5) / E^5 * b * r / rss^2 * (Tc * cos(TA) + Ts * sin(TA));
H6 = 45 * RE^4 * J2^2 * mu * a^2 / r^2 * (6*e^2*(15*s^2 - 14)*s^2*cos(2*g) - 3 * e^2 * (5*s^4+8*s^2-8) + 84*s^4-168*s^2+80) / 8 / a^5 / (e^2-1)^3 / E^6;

if getHtilde == 1
    Htilde = TotalHamiltonian([H0, 0, H2, H3, 0, H5, H6],E);
    Htilde = substituteFunctions(Htilde, rfun, TAfun, Eccfun);
    for j = 1:8
        if j < 5
            FullHDynamics(j,1) = diff(Htilde,Xx(j));
        else
            FullHDynamics(j,1) = -diff(Htilde,Xx(j));
        end
    end
    FullHDynamics = [FullHDynamics(1:3);FullHDynamics(5:7);FullHDynamics(4);FullHDynamics(8)];
    FullHDynamics = substitutePartials(FullHDynamics, subPartialsInputs);
    if writeNewFile == 1
        FullHDynamics = matlabFunction(FullHDynamics,'File','HtildeDynamicsFunction');
    end
end

%% Transformed Hamiltonians from Inkwan Park

K00 = -mu / 2 / a;
K02 = 2 * Ks * nu / E^2;
K03 = 3 * RE^2 * J2 * mu * (3 * s^2 - 2) / 2 / a^3 / eta^3 / E^3;
K05 = - 180 * a * b * e * Tc / E^5;
K06 = - 45 * RE^4 * J2^2 * mu  / 8 / a^5 / eta^9 / E^6 * (6 * e^2 * eta^2 * s^2 * (15 * s^2 - 14) * cos(2*g) ...
    + e^2 * (- (3 * eta^2 * (5 * s^4 + 8 * s^2 - 8) + 4 * (2-3*s^2)^2)) ...
    + 4 * (3 * eta^3 * (2-3*s^2)^2 + eta^2 * (21 * s^4 - 42 * s^2 + 20) + (2-3*s^2)^2)); 
% K06 = 0;
if getK == 1
    Kn = TotalHamiltonian([K00, 0, K02, K03, 0, K05, K06],E);
    Kn = substituteFunctions(Kn, rfun, TAfun, Eccfun);
    for j = 1:8
        if j < 5
            MeanDynamics(j,1) = diff(Kn,Xx(j));
        else
            MeanDynamics(j,1) = -diff(Kn,Xx(j));
        end
    end
    
    MeanDynamics = substitutePartials(MeanDynamics, subPartialsInputs);
    MeanDynamics = [MeanDynamics(1:3);MeanDynamics(5:7);MeanDynamics(4);MeanDynamics(8)];
    if writeNewFile == 1
        MeanDynamics = matlabFunction(MeanDynamics,'File','MeanDynamicsFunction');
    end

end

%% Normalizing generating functions from Inkwan Park
W3n = -factorial(3) / E^3 * L^3 * RE^2 * J2 / 4 / mu / a^3 / eta^3 * ((3*s^2-2) * (TA + e*sin(TA) - l) + ...
        (-3*s^2)*(1/2*sin(2*(TA+g))+1/2*e*sin(TA+2*g) + 1/6*e*sin(3*TA+2*g)));
W5n = -30 * a * b * (4 * e^2 * Tc * sin(Ecc) + e * (eta * Ts * cos(2*Ecc) - Tc * (6*Ecc+sin(2*Ecc)-6*l)) - 4 * eta * Ts * cos(Ecc) + 4 * Tc * sin(Ecc)) / n / E^5;
W6n = -45 * RE^4 * J2^2 * n / 8 / a^2 / e / eta^9 / E^6 * ...
        (-6 * e^3 * eta^2 * s^2 * (15*s^2 - 14) * (TA-l) * cos(2*g) ...
        + 8 * (2-3*s^2)^2*(e^2 + eta^3 - 1) * sin(TA) ...
        + e * ( 2 * (2-3*s^2)^2 * (e^2 + eta^3 - 1) * sin(2*TA) ...
        + (TA-l) * (e^2 * (3*eta^2*(5*s^4 + 8*s^2 - 8) + 4*(2-3*s^2)^2) ...
        - 4 * (eta^2 * (21*s^4 - 42 * s^2 + 20) + (2-3*s^2)^2))));
if getKLW1 == 1
    W3n = substituteFunctions(W3n, rfun, TAfun, Eccfun);
    W5n = substituteFunctions(W5n, rfun, TAfun, Eccfun);
    W6n = substituteFunctions(W6n, rfun, TAfun, Eccfun);
    H0s = substituteFunctions(H0, rfun, TAfun, Eccfun);
    
    % Find variable terms
    LoW3 = - diff(H0s,L) * diff(W3n,l);
    LoW5 = - diff(H0s,L) * diff(W5n,l);
    LoW6 = - diff(H0s,L) * diff(W6n,l);
    
    LoW3 = substitutePartials(LoW3, subPartialsInputs);
    LoW5 = substitutePartials(LoW5, subPartialsInputs);
    LoW6 = substitutePartials(LoW6, subPartialsInputs);
      
    LW = TotalHamiltonian([0,0,0,LoW3,0,LoW5,LoW6],E);
    LW = substituteFunctions(LW, rfun, TAfun, Eccfun);
    
    % Mean Dynamics 
    for j = 1:8
        if j < 5
            ShortPeriodDynamics(j,1) = diff(LW,Xx(j));
        else
            ShortPeriodDynamics(j,1) = -diff(LW,Xx(j));
        end
    end
    ShortPeriodDynamics = substitutePartials(ShortPeriodDynamics, subPartialsInputs);
    ShortPeriodDynamics = [ShortPeriodDynamics(1:3);ShortPeriodDynamics(5:7);ShortPeriodDynamics(4);ShortPeriodDynamics(8)];
    
    % Full dynamics
    Ht03 = K03 + LoW3;
    Ht05 = K05 + LoW5;
    Ht06 = K06 + LoW6;
  
    Knt = TotalHamiltonian([K00,0,K02,Ht03,0,Ht05,Ht06],E);
    Knt = substituteFunctions(Knt, rfun, TAfun, Eccfun);
    
    for j = 1:8
        if j < 5
            FullDynamics(j,1) = diff(Knt,Xx(j));
        else
            FullDynamics(j,1) = -diff(Knt,Xx(j));
        end
    end
    
    FullDynamics = substitutePartials(FullDynamics, subPartialsInputs);
    FullDynamics = [FullDynamics(1:3);FullDynamics(5:7);FullDynamics(4);FullDynamics(8)];
    if writeNewFile == 1
        ShortPeriodDynamics = matlabFunction(ShortPeriodDynamics,'File','ShortPeriodDynamicsFunction');
        FullDynamics = matlabFunction(FullDynamics,'File','FullDynamicsFunction');
    end

end

if getKLW2 == 1
    SP0 = H0 - K00;
    SP2 = H2 - K02;
    SP3 = H3 - K03;
    SP5 = H5 - K05;
    SP6 = LW6;
    SP = TotalHamiltonian([SP0,0,SP2,SP3,0,SP5,SP6],E);
    SPfun = substituteFunctions(SP, rfun, TAfun, Eccfun);
    
    % Mean Dynamics 
    for j = 1:8
        if j < 5
            ShortPeriodDynamics(j,1) = diff(SPfun,Xx(j));
        else
            ShortPeriodDynamics(j,1) = -diff(SPfun,Xx(j));
        end
    end
    ShortPeriodDynamics = substitutePartials(ShortPeriodDynamics, subPartialsInputs);
    ShortPeriodDynamics = [ShortPeriodDynamics(1:3);ShortPeriodDynamics(5:7);ShortPeriodDynamics(4);ShortPeriodDynamics(8)]
    if writeNewFile == 1
        ShortPeriodDynamics = matlabFunction(ShortPeriodDynamics,'File','ShortPeriodDynamicsFunction');    
    end
end

%%

if flagGetInitialDelOffset == 1
    W3n = substituteFunctions(W3n, rfun, TAfun, Eccfun);
    W5n = substituteFunctions(W5n, rfun, TAfun, Eccfun);
    W6n = substituteFunctions(W6n, rfun, TAfun, Eccfun);
    order = 6;
    q = [l;g;h;ks];
    Q = [L;G;H;Ks];
    W = [0;0;W3n;0;W5n;W6n];
    for N = 1:order
        for i = 1:4
            dWndYNi = diff(W(N),Q(i));
            dWndyNi = diff(W(N),q(i));
            dWndY(N,i) = substitutePartials(dWndYNi, subPartialsInputs);
            dWndy(N,i) = substitutePartials(dWndyNi, subPartialsInputs);
        end
    end
    
    % Get yon and Yon
    for i = 1:4
        SUM1 = 0;
        SUM2 = 0;
        clearvars yo Yo GjYoNmj GjyoNmj CLyoGyo CLYoGYo LjyoNmj LjYoNmj xon yoN YoN Xon
        for N = 1:order
            if N == 1
                CGy1 = 0;
                CGY1 = 0;
                CGy2 = 0;
                CGY2 = 0;
            else
                sum1 = 0;
                sum2 = 0;
                sum3 = 0;
                sum4 = 0;
                for j = 1:N-1
                    if j == 1 || j == 2
                        CLyoGyo = 0;
                        CLYoGYo = 0;
                    else
                        sum01 = 0;
                        sum02 = 0;
                        for m = 0:j-2
                            LGyo = diff(Gyo(j-m-1,N-j),q(i)) * diff(W(m+1),Q(i)) - diff(Gyo(j-m-1,N-j),Q(i)) * diff(W(m+1),q(i));
                            LGYo = diff(GYo(j-m-1,N-j),q(i)) * diff(W(m+1),Q(i)) - diff(GYo(j-m-1,N-j),Q(i)) * diff(W(m+1),q(i));
                            LGyo = substitutePartials(LGyo, subPartialsInputs);
                            LGYo = substitutePartials(LGYo, subPartialsInputs);
                            sum01 = sum01 + nchoosek(j-1,m) * LGyo;
                            sum02 = sum02 + nchoosek(j-1,m) * LGYo;
                        end
                        CLyoGyo = sum01;
                        CLYoGYo = sum02;
                    end
                    LjyoNmj = diff(yo(N-j),q(i)) * diff(W(j),Q(i)) - diff(yo(N-j),Q(i)) * diff(W(j),q(i));
                    LjYoNmj = diff(Yo(N-j),q(i)) * diff(W(j),Q(i)) - diff(Yo(N-j),Q(i)) * diff(W(j),q(i));
                    LjyoNmj = substitutePartials(LjyoNmj, subPartialsInputs);
                    LjYoNmj = substitutePartials(LjYoNmj, subPartialsInputs);
                    GjyoNmj = LjyoNmj - CLyoGyo;
                    GjYoNmj = LjYoNmj - CLYoGYo;
                    Gyo(j,N-j) = substituteFunctions(GjyoNmj, rfun, TAfun, Eccfun);
                    GYo(j,N-j) = substituteFunctions(GjYoNmj, rfun, TAfun, Eccfun);
                    sum1 = sum1 + nchoosek(N-1,j) * GjyoNmj;
                    sum2 = sum2 + nchoosek(N-1,j) * GjYoNmj;
                    sum3 = sum3 + nchoosek(N,j) * GjyoNmj;
                    sum4 = sum4 + nchoosek(N,j) * GjYoNmj;
                end
                CGy1 = sum1;
                CGY1 = sum2;
                CGy2 = sum3;
                CGY2 = sum4;
            end
            
            yoN = dWndY(N,i) + CGy1;
            YoN = -dWndy(N,i) + CGY1;
            xon = -yoN + CGy2;
            Xon = -YoN + CGY2;
            yo(N) = substituteFunctions(yoN, rfun, TAfun, Eccfun);
            Yo(N) = substituteFunctions(YoN, rfun, TAfun, Eccfun);
            SUM1 = SUM1 + E^N / factorial(N) * xon;
            SUM2 = SUM2 + E^N / factorial(N) * Xon;
        end
        yOffset(i) = SUM1;
        YOffset(i) = SUM2;
    
    end
    
    getInitialDelOffset = [yOffset,YOffset]';
    getInitialDelOffset = [getInitialDelOffset(1:3);getInitialDelOffset(5:7);getInitialDelOffset(4);getInitialDelOffset(8)];
    if writeNewFile == 1
        getInitialDelOffset = matlabFunction(getInitialDelOffset,'File','getInitialDelOffset');
    end
end

if flagGetFinalDelOffset == 1
    W3n = substituteFunctions(W3n, rfun, TAfun, Eccfun);
    W5n = substituteFunctions(W5n, rfun, TAfun, Eccfun);
    W6n = substituteFunctions(W6n, rfun, TAfun, Eccfun);
    order = 6;
    q = [l;g;h;ks];
    Q = [L;G;H;Ks];
    W = [0;0;W3n;0;W5n;W6n];
    for N = 1:order
        for i = 1:4
            dWndYNi = diff(W(N),Q(i));
            dWndyNi = diff(W(N),q(i));
            dWndY(N,i) = substitutePartials(dWndYNi, subPartialsInputs);
            dWndy(N,i) = substitutePartials(dWndyNi, subPartialsInputs);
        end
    end
    
    % Get yon and Yon
    for i = 1:4
        SUM1 = 0;
        SUM2 = 0;
        clearvars yo Yo GjYoNmj GjyoNmj CLyoGyo CLYoGYo LjyoNmj LjYoNmj xon yoN YoN Xon
        for N = 1:order
            if N == 1
                CGy1 = 0;
                CGY1 = 0;
                CGy2 = 0;
                CGY2 = 0;
            else
                sum1 = 0;
                sum2 = 0;
                for j = 1:N-1
                    if j == 1 || j == 2
                        CLyoGyo = 0;
                        CLYoGYo = 0;
                    else
                        sum01 = 0;
                        sum02 = 0;
                        for m = 0:j-2
                            LGyo = diff(Gyo(j-m-1,N-j),q(i)) * diff(W(m+1),Q(i)) - diff(Gyo(j-m-1,N-j),Q(i)) * diff(W(m+1),q(i));
                            LGYo = diff(GYo(j-m-1,N-j),q(i)) * diff(W(m+1),Q(i)) - diff(GYo(j-m-1,N-j),Q(i)) * diff(W(m+1),q(i));
                            LGyo = substitutePartials(LGyo, subPartialsInputs);
                            LGYo = substitutePartials(LGYo, subPartialsInputs);
                            sum01 = sum01 + nchoosek(j-1,m) * LGyo;
                            sum02 = sum02 + nchoosek(j-1,m) * LGYo;
                        end
                        CLyoGyo = sum01;
                        CLYoGYo = sum02;
                    end
                    LjyoNmj = diff(yo(N-j),q(i)) * diff(W(j),Q(i)) - diff(yo(N-j),Q(i)) * diff(W(j),q(i));
                    LjYoNmj = diff(Yo(N-j),q(i)) * diff(W(j),Q(i)) - diff(Yo(N-j),Q(i)) * diff(W(j),q(i));
                    LjyoNmj = substitutePartials(LjyoNmj, subPartialsInputs);
                    LjYoNmj = substitutePartials(LjYoNmj, subPartialsInputs);
                    GjyoNmj = LjyoNmj - CLyoGyo;
                    GjYoNmj = LjYoNmj - CLYoGYo;
                    Gyo(j,N-j) = substituteFunctions(GjyoNmj, rfun, TAfun, Eccfun);
                    GYo(j,N-j) = substituteFunctions(GjYoNmj, rfun, TAfun, Eccfun);
                    sum1 = sum1 + nchoosek(N-1,j) * GjyoNmj;
                    sum2 = sum2 + nchoosek(N-1,j) * GjYoNmj;
                end
                CGy1 = sum1;
                CGY1 = sum2;
            end
            
            yoN = dWndY(N,i) + CGy1;
            YoN = -dWndy(N,i) + CGY1;
            SUM1 = SUM1 + E^N / factorial(N) * yoN;
            SUM2 = SUM2 + E^N / factorial(N) * YoN;
            yo(N) = substituteFunctions(yoN, rfun, TAfun, Eccfun);
            Yo(N) = substituteFunctions(YoN, rfun, TAfun, Eccfun);
        end
        xOffset(i) = SUM1;
        XOffset(i) = SUM2;
    
    end
    
    getFinalDelOffset = [xOffset,XOffset]';
    getFinalDelOffset = [getFinalDelOffset(1:3);getFinalDelOffset(5:7);getFinalDelOffset(4);getFinalDelOffset(8)];
    if writeNewFile == 1
        FinalDelOffset = matlabFunction(getFinalDelOffset,'File','getFinalDelOffset');
    end
end

%% Functions
function OUT = TotalHamiltonian(Hs,E)
    Htot = 0;
    
    for j = 1:length(Hs)
        if j == 1
            Htoti(j) = Hs(1);
        else
            Htoti(j) = E^(j-1) / factorial(j-1) * Hs(j);
        end
        Htot = Htot + Htoti(j);
    end
    OUT = Htot;
end


