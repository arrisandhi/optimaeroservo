%% DEFINE DESIGN VECTOR
% Design vector contains three coefficient[a;b;c] of quadratic mass
% distribution  RA(i) = ax(i)^2+bx(i)+c
tic
X0 = [0;0;0;0;35.71] ;
X0 = [-0.00115817543250154;0.0519394788147191;-0.579443662404584;0.852629693983976;35.2871701074614];

%% SIMULATE GUST RESPONSE OF ORIGINAL DESIGN
main_pre_MF05
main_processing
[u,x,y] = main_LQR(ssLinRed,StateWORBRed,SimDef,GustDef);
zDispWORBRed = calc_disprot(u,x,y,BeamDef,SimDef,trimstate,StateWORBRed,StateWORB,Nmodes,BeamSeg,Mode,Ti,RVnum);

% Calculate total wing mass
BaseRA = max(abs(BeamProp.RA));
% Calculate strain energy of initial configuration
maxMomOri = max(abs(y(1,:)));


%% OPTIMISATION
[history,searchdir] = runfmincon_MF05_TM_MMax(maxMomOri,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,X0,BaseRA);
timeSimulation = toc;
filename = 'MF05_TM_10MMax.mat';
save(filename);
