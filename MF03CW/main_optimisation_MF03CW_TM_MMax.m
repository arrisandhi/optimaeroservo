%% DEFINE DESIGN VECTOR
% Design vector contains three coefficient of quadratic mass
% distribution i.e. RA(i) = ax(i)^2+bx(i)+c and three coefficient of
% control weight Q and R. X = [a;b;c;Q;R]
tic
X0 = [0;0;35.71;20;50] ;

%% SIMULATE GUST RESPONSE OF ORIGINAL DESIGN
main_pre_MF03CW
main_processing
CW = X0(4:5);       % Control weight
[u,x,y] = main_LQR_MF03CW(ssLinRed,StateWORBRed,SimDef,GustDef,CW);
zDispWORBRed = calc_disprot(u,x,y,BeamDef,SimDef,trimstate,StateWORBRed,StateWORB,Nmodes,BeamSeg,Mode,Ti,RVnum);

% Calculate total wing mass
BaseRA = max(abs(BeamProp.RA));
% Calculate strain energy of initial configuration
maxMomOri = max(abs(y(1,:)));


%% OPTIMISATION
% X0 = [-0.156234058002392;5.04130330520589;10;20;50] ;
% X0 = [-0.0627004898652805;1.82281229433182;35.5928140469366;20;50];
[history,searchdir] = runfmincon_MF03CW_TM_MMax(maxMomOri,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,X0,BaseRA);
timeSimulation = toc;
filename = 'MF03CW_TM_10MMax_X0Goland.mat';
save(filename);
