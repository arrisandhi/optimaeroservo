function BeamProp = main_beamprop(BeamDef)
%% THIS CODE RETURNS THE PROPERTIES OF WING
% INPUTS:
% BeamDef:
%   h1      & 1x1 & Height of idealized wing section
%   h2      & 1x1 & Width of idealized wing section
%   rho     & 1x1 & Density of wing 
%   Emod    & 1x1 &
%   Gmod    & 1x1 &
%   Ndivs   % 1x1 % Number of wing element

h1 = BeamDef.h1;
h2 = BeamDef.h2;
rho = BeamDef.rho;
Emod = BeamDef.Emod;
Gmod = BeamDef.Gmod;
Ndivs = BeamDef.Ndivs;

%% FOR BEAM PROPERTIES ACCORDING TO PROPERTIES STATED IN MAIN_PRE
% RA      = rho*h1*h2;                        % Mass per unit length 
% RI1     = rho*h1*h2*(h1*h1+h2*h2)/12;       % Sectional moment inertia per unit length in x direction
% RI2     = rho*h1*h2^3/12;                   % Sectional moment inertia per length unit in y direction
% RI3     = rho*h1^3*h2/12;                   % Sectional moment nertia per length unit in z direction
% EA      = Emod*h1*h2;                       % Young Modulus 
% EI1     = h1^3*h2/12*Emod;                  % Torsional modulus in z direction
% EI2     = h1*h2^3/12*Emod;                  % Torsional modulus in y direction 
% GJ      = 0.229*h1*h2*h1*h2*Gmod;           % Torsional modulus in x direction

%% FOR USING GOLAND WING. USE THE PROPERTIES BELOW
RA = 35.71*ones(Ndivs,1);
RI1 = 8.64*ones(Ndivs,1);
RI2 = RI1;
RI3 = RI1;
EA  = 1e10*ones(Ndivs,1);
GJ  = 0.99e6*ones(Ndivs,1);
EI2 = 9.77e6*ones(Ndivs,1);
EI1 = 1e8*ones(Ndivs,1);
%Note that for Goland wing the cg is 43% chord, therefore we should add this following code in beamconvergence_Arri Code 
% StandardM(3,4) = -0.2*0.9144*RA;
% StandardM(4,3) = StandardM(3,4);



%% Baseline configuration
RA  = [8.93*ones(19,1);5]; 
RI1 = 4.15*ones(Ndivs,1);
RI2 = 0.69*ones(Ndivs,1);
RI3 = 3.46*ones(Ndivs,1);
EA  = 1e10*ones(Ndivs,1);
GJ  = 1.65e5*ones(Ndivs,1);
EI2 = 1.03e6*ones(Ndivs,1);
EI1 = 1.24e7*ones(Ndivs,1);


% RA  = 5*ones(Ndivs,1); 
% RA = [7.34E+00;9.14E+00;7.04E+00;5.11E+00;6.61E+00;5.01E+00;9.63E+00;8.85E+00;8.64E+00;5.00E+00;5.06E+00;5.17E+00;1.31E+01;1.55E+01;5.91E+00;8.06E+00;1.02E+01;1.17E+01;1.25E+01;8.53E+00];

% Optimised config, Flap 1-6
% RA = [8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.94202510527019;8.28973247091385;18.1165989860157];

% Optimised config, Flap 1-11
% RA = [9.00796898820397;9.00796898820435;9.00796898820409;9.00796898820428;9.00796898820435;9.00796898820430;9.00796898820425;9.00796898820447;9.00796898820432;9.00796898820416;9.00796898820401;9.00796898820413;9.00796898820403;9.00796898820401;9.00796898820415;9.00796898820419;9.00796898820406;9.00796898820403;6.63459349487089;7.67218420719878];

% Optimised config, Flap 1-21
% RA = [5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5.00115740643061;8.28628472040847];

% Optimised config, Flap 6-11
% RA = [9.37076290336696;9.37076290336506;9.37076290336665;9.37076290513164;9.37076290336665;9.37076290513163;9.37076290513323;9.37076290336505;9.37076290336664;9.37076290336505;9.37076290336664;9.37076290336505;9.37076290336664;9.37076290336504;9.37076290336664;9.37076290336504;9.37076290336664;9.37076290336504;25.7757732753082;7.99581231319768];

% Optimised config, Flap 11-16
% RA = [5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;8.97040328083288;11.7947945483907];

% Optimised config, Flap 11-21
% RA = [5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5.00023659774666;18.5656403903854];

% Optimised config, Flap 16-21
% RA = [5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;7.24096100295248;7.59198502779831];
% RA = [5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5;5.71436143029214;12.0221593881705];

%% ASSIGN TO BEAMPROP DATA STRUCTURE
BeamProp.RA  = RA;
BeamProp.RI1 = RI1;
BeamProp.RI2 = RI2;
BeamProp.RI3 = RI3;
BeamProp.EA  = EA;
BeamProp.EI1 = EI1;
BeamProp.EI2 = EI2;
BeamProp.GJ  = GJ;
end
