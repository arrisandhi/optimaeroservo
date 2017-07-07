function [Cost] = main_optimCost_MF02_TM_MMax(X0,BeamDef,BaseRA)
% INPUTS:
    % X0     & 2x1          & Design vector
    % NodeL  & NumNode x 1  & x-coordinate of nodes
    % Ndivs  & 1x1          & Number of wing elements 
    % BaseRA & 1x1          & Maximum mass density of initial wing configuration

% OUTPUTS:
    % Cost   & 1x1          & Optimisation cost

a = X0(1);
b = X0(2);
NodeL = BeamDef.NodeL;
Ndivs = BeamDef.Ndivs;

RA = zeros(Ndivs,1);
for i = 1:Ndivs;
    RA(i) = a*NodeL(i)+b;
end

%% COST CALCULATION
% Normalised design varible value by max mass density bound

NormRA = RA/BaseRA;
Cost = sum(NormRA);
end