%% DEFINE DESIGN VECTOR
% Design vector contains three coefficient[a;b;c] of quadratic mass
% distribution  RA(i) = ax(i)^2+bx(i)+c
tic
X0 = [0;0;35.71] ;
% X0 = [-1.562341e-01;5.041303e+00;8.271621e+00];
%% SIMULATE GUST RESPONSE OF ORIGINAL DESIGN
main_pre_MF03_TM_MMMS
main_processing
[u,x,y] = main_LQR(ssLinRed,StateWORBRed,SimDef,GustDef);
zDispWORBRed = calc_disprot(u,x,y,BeamDef,SimDef,trimstate,StateWORBRed,StateWORB,Nmodes,BeamSeg,Mode,Ti,RVnum);

% Calculate total wing mass
BaseRA = max(abs(BeamProp.RA));
% Calculate strain energy of initial configuration
maxMomOri = max(abs(y(1,:)));


%% OPTIMISATION
X0 = [0;-2.06;64.64936];
[history,searchdir] = runfmincon_MF03_TM_MMMS(maxMomOri,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,X0,BaseRA);
timeSimulation = toc;
filename = 'MF03_TM_10MMMS.mat';
save(filename);
