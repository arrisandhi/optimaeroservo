function BeamProp = main_beamprop_MF02(BeamDef,ShapeMass)
%% THIS CODE RETURNS THE PROPERTIES OF WING
% INPUTS:
    % ShapeMass   & 2x1 & Design vector
    % BeamDef:
        % NodeL       & 1x1 & Length of wing element 
        % Ndivs       % 1x1 % Number of wing element
% OUTPUTS:
    % BeamProp    & struct & Beam properties  
        % RA        : Mass per unit length
        % RI1       : Sectional moment inertia per unit length(x direction)
        % RI2       : Sectional moment inertia per unit length(y direction)
        % RI3       : Sectional moment inertia per unit length(z direction)
        % EA        : Young Modulus in x,y,z direction
        % EI1       : In-plane bending stiffness
        % EI2       : Out-plane bending stiffness 
        % GJ        : Twisting stiffness
        
a = ShapeMass(1);       % Gradient of mass density distribution
b = ShapeMass(2);       % Coefficient 
NodeL = BeamDef.NodeL;  
Ndivs = BeamDef.Ndivs;

RA = zeros(Ndivs,1);
for i = 1:Ndivs;
    RA(i) = a*NodeL(i)+b;
end
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

% Store properties in BeamProp structure array.
BeamProp.RA  = RA;
BeamProp.RI1 = RI1;
BeamProp.RI2 = RI2;
BeamProp.RI3 = RI3;
BeamProp.EA  = EA;
BeamProp.EI1 = EI1;
BeamProp.EI2 = EI2;
BeamProp.GJ  = GJ;
end
