function [Cost] = main_optimCost_MF04_TM_MMax(X0,BeamDef,BaseRA)
a = X0(1);
b = X0(2);
c = X0(3);
d = X0(4);
NodeL = BeamDef.NodeL;
Ndivs = BeamDef.Ndivs;

RA = zeros(Ndivs,1);
for i = 1:Ndivs;
    RA(i) = a*NodeL(i)^3+b*NodeL(i).^2+c*NodeL(i)+d;
end

%% COST CALCULATION
% Normalised design varible value by max mass density bound

NormRA = RA/BaseRA;

Cost = sum(NormRA);
end