function [cnlin, cnleq] = confun_MF02_TM_MMax(X0,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,maxMomOri)
a = X0(1);
b = X0(2);
NodeL = BeamDef.NodeL;
Ndivs = BeamDef.Ndivs;

RA = zeros(Ndivs,1);
for i = 1:Ndivs;
    RA(i) = a*NodeL(i)+b;
end
BeamProp.RA = RA;


% RERUN main_processing FOR NEW WING CONFIGURATION
%% OBTAIN K AND M MATRIX
[StructMat] = beamconvergence(BeamDef,BeamProp);

%% CALCULATE MATRIX COEFFICIENT
[MC,~,AeDef,Mode,btoV] = guyancomp(Nmodes,BeamSeg,BeamDef,AeDef,StructMat);

%% ASSEMBLY TRIMSTATE
trimstate = assy_trimstate(AeDef,Nmodes,btoV);

%% ASSEMBLY EIGMAT(NOT INCLUDING WD MATRIX AND HAS NOT BEEN MULTIPLIED BY AMAT_INVERSE YET)
EigMatFull = systemdxeigconvergence(MC,AeDef,trimstate,btoV,Nmodes);

%% CALCULATE INPUT AND OUTPUT MATRIX
[SBG,SCG] = systemdxioGust(MC,BeamDef,GustDef,AeDef,BeamSeg,Mode,trimstate,Nmodes);

%% FOR CANTILEVER WING, STATES RELATED TO RIGID BODY WING MUST BE EXCLUDED USING THE FOLLOWING
% LINE. WORB(WITHOUT RIGID BODY)
[wingSSG,StateWORB] = excludeRB(MC,EigMatFull,SBG,SCG,trimstate,RVnum,Nmodes);

%% FOR FREE FLYING WING, IT IS NOT REQUIRED TO EXCLUDE RIGID BODY STATE.
% USE THIS LINE TO OBTAIN DYNAMIC MATRIX AND COMMENT THE PREVIOUS CODE SECTION
% EigMatFull(1:2*neig+6,:) = AMat\EigMatFull(1:2*neig+6,:);
% wingSSG = ss(EigMatFull,SBG,SCG,zeros(size(SCGWORB,1),size(SBGWORB,2));

%% MODEL ORDER REDUCTION
[ssLinRed,StateWORBRed,Ti] = MOR(wingSSG,StateWORB);

%% SIMULATION
fprintf(1, 'Run LQR.\n');
% [u,du,x,y] = main_MPC(ssLinRed,StateWORBRed,SimDef,GustDef);
[u,x,y,Klqr] = main_LQR(ssLinRed,StateWORBRed,SimDef,GustDef);
maxMomRoot = max(abs(y(1,:)));

%% COST CALCULATION
% Nonlinear inequality constraints. 
% In this case we want the optimised wing has 5% MaxMom reduction
cMom = maxMomRoot - 0.95*maxMomOri;

% Stability constraint
eigFlutter = max(real(eig(ssLinRed.a)));
if eigFlutter > 0
    cFS = 1;
else
    cFS = 0;
end

cnlin = [cMom ; cFS];

% Nonlinear equality constraints
cnleq = [];
end