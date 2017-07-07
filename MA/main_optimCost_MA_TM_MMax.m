function [Cost] = main_optimCost_MA_TM_MMax(X0,BeamDef,BaseRA)
RA = X0;

%% COST CALCULATION
% Normalised design varible value by max mass density bound

NormRA = RA/BaseRA;

Cost = sum(NormRA);
end