function BeamProp = main_beamprop_MA(BeamDef,X0)
%% THIS CODE RETURNS THE PROPERTIES OF WING
% INPUTS:
% BeamDef:
%   h1      & 1x1 & Height of idealized wing section
%   h2      & 1x1 & Width of idealized wing section
%   rho     & 1x1 & Density of wing 
%   Emod    & 1x1 &
%   Gmod    & 1x1 &
%   Ndivs   % 1x1 % Number of wing element
Ndivs = BeamDef.Ndivs;

RA = X0;
RI1 = 8.64*ones(Ndivs,1);
RI2 = RI1;
RI3 = RI1;
EA  = 1e10*ones(Ndivs,1);
GJ  = 0.99e6*ones(Ndivs,1);
EI2 = 9.77e6*ones(Ndivs,1);
EI1 = 1e8*ones(Ndivs,1);


%Note that for Modified Goland wing the cg is 23%, therefore we should add this following code in beamconvergence_Arri Code 
% StandardM(3,4) = 0.2*0.9144*RA;
% StandardM(4,3) = StandardM(3,4);

BeamProp.RA  = RA;
BeamProp.RI1 = RI1;
BeamProp.RI2 = RI2;
BeamProp.RI3 = RI3;
BeamProp.EA  = EA;
BeamProp.EI1 = EI1;
BeamProp.EI2 = EI2;
BeamProp.GJ  = GJ;
end
