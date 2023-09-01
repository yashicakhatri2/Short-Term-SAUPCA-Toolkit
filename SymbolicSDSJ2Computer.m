% Symbolic Jacobian Computer SDS
function OUT = SymbolicSDSJ2Computer(mu, RE, J2, order,dyn,writeNewFile)

syms l g h L G H ks Ks hs cs sSun MU
syms J2 b mu RE nu r rfun(L,G,H,l) crLpart crGpart crHpart crlpart TAfun(L,G, H,l) TA cfLpart cfGpart cflpart cfHpart Ecc Eccfun(L,G,H,l) cELpart cEGpart cElpart cEHpart E

x = [l;g;h;ks]; X = [L;G;H;Ks];

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

crlpartVal = drde * dedl + drdM * dMdl; % + drda * dadl;
crLpartVal = drde * dedL + drdM * dMdL; % + drda * dadL;
crGpartVal = drde * dedG + drdM * dMdG; % + drda * dadG;
crHpartVal = drde * dedH + drdM * dMdH; % + drda * dadH;

cElpartVal = dEde * dedl + dEdM * dMdl;
cELpartVal = dEde * dedL + dEdM * dMdL;
cEGpartVal = dEde * dedG + dEdM * dMdG;
cEHpartVal = dEde * dedH + dEdM * dMdH;

cflpartVal = dfde * dedl + dfdM * dMdl;
cfLpartVal = dfde * dedL + dfdM * dMdL; 
cfGpartVal = dfde * dedG + dfdM * dMdG; 
cfHpartVal = dfde * dedH + dfdM * dMdH; 

subPartialsInputs = [r, crLpart, crGpart, crHpart, crlpart, cfLpart, cfGpart, cfHpart, cflpart, TA, cELpart, cEGpart, cElpart, cEHpart, Ecc,crLpartVal,crGpartVal,crHpartVal,crlpartVal,cfLpartVal,cfGpartVal,cfHpartVal, cflpartVal,cELpartVal,cEGpartVal,cElpartVal,cEHpartVal];

if dyn == "mean"
    DOT = MeanDynamicsFunction(G,H,J2,L,RE,b,cs,g,h,hs,ks,mu,nu,sSun);
elseif dyn == "SP"
    DOT = ShortPeriodDynamicsFunction(Ecc,G,H,J2,L,RE,TA,b,cs,g,h,hs,ks,mu,r,sSun);
end
DOT = substituteFunctions(DOT, rfun, TAfun, Eccfun);

State = [l;g;h;L;G;H;ks;Ks];

%% First order partials calcs
J = jacobian(DOT, State);
J = substitutePartials(J, subPartialsInputs);

if order == 1
    OUT = reshape(J',64,1);
elseif order == 2
    Jsub = substituteFunctions(J, rfun, TAfun, Eccfun);
    dim = 8;
    for i = 1:dim
        for j = 1:dim
            for n = 1:dim
                ANew(i,j,n) = diff(Jsub(i,j),State(n));
            end
        end
    end
    ANew = substitutePartials(ANew, subPartialsInputs);
    for n = 1:dim
        Amat(:,:) = ANew(:,:,n);
        ANewCol(:,n) = reshape(Amat',dim^2,1);
    end
    AFinal = reshape(ANewCol,dim^3,1);
    OUT = AFinal;
end
if writeNewFile == 1
    getJOrder2 = matlabFunction(OUT,'File','getJOrder2');
end

end
