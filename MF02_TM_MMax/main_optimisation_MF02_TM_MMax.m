%% DEFINE DESIGN VECTOR
% Design vector contains two coefficient(a and b) of linear mass
% distribution  RA(i) = ax(i)+b
tic
X0 = [0;35.71] ;

%% SIMULATE GUST RESPONSE OF ORIGINAL DESIGN
main_pre_MF02
main_processing
[u,x,y,Klqr] = main_LQR(ssLinRed,StateWORBRed,SimDef,GustDef);
% [u,du,x,y] = main_MPC(ssLinRed,StateWORBRed,SimDef,GustDef); % uncomment this for using MPC

% Calculate maximum mass density of initial configuration for cost optimisation calculation
BaseRA = max(abs(BeamProp.RA));
% Calculate maximum moment of initial configuration for optimisation reference
maxMomOri = max(abs(y(1,:)));

%% OPTIMISATION
[history,searchdir] = runfmincon_MF02_TM_MMax(maxMomOri,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,ShapeMass,BaseRA);
timeSimulation = toc;
filename = 'MF02_TM_05MMax.mat';
save(filename);
