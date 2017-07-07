%% DEFINE DESIGN VECTOR
% Design vector contains all mass per unit length for each wing elements
tic
% X0 = 35.71*ones(20,1) ;
% 10% MMax Reduction
X0 = [5.01675734453503;22.2853254823265;18.5035329154662;33.0614376878825;33.7895450719233;31.6343720289448;31.1402250525112;32.7877447464512;40.4947004584155;40.8958938651335;44.8841975552928;47.0647969499637;51.1645075087336;51.1413933215986;56.3035675228940;48.1941669951215;40.0561452269584;33.4256856989902;28.6594109422595;25.3646667232919];

%% SIMULATE GUST RESPONSE OF ORIGINAL DESING
main_pre_MA
main_processing
[u,x,y] = main_LQR(ssLinRed,StateWORBRed,SimDef,GustDef);
% Calculate total wing mass
BaseRA = max(abs(BeamProp.RA));
% Calculate strain energy of initial configuration
maxMomOri = max(abs(y(1,:)));


%% OPTIMISATION
[history,searchdir] = runfmincon_MA_TM_MMax(maxMomOri,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,X0,BaseRA);
timeSimulation = toc;
filename = 'MA_TM_10MMax.mat';
save(filename);
