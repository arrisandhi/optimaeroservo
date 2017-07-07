function [Cost] = main_optimCost_MF05_TM_MMax(X0,BeamDef,BaseRA)
a = X0(1);
b = X0(2);
c = X0(3);
d = X0(4);
e = X0(5);
NodeL = BeamDef.NodeL;
Ndivs = BeamDef.Ndivs;

RA = zeros(Ndivs,1);
for i = 1:Ndivs;
    RA(i) = a*NodeL(i)^4+b*NodeL(i).^3+c*NodeL(i).^2+d*NodeL(i)+e;
end

%% COST CALCULATION
% Normalised design varible value by max mass density bound

NormRA = RA/BaseRA;

Cost = sum(NormRA);
end