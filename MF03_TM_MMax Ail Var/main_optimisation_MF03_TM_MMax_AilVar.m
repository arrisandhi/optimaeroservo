%% DEFINE DESIGN VECTOR
% Design vector contains three coefficient[a;b;c] of quadratic mass
% distribution  RA(i) = ax(i)^2+bx(i)+c
tic
X0 = [0;0;35.71] ;
% X0 = [0.0333492240548254;0.118357385527683;32.0613494339622];   %10Ail_1
% TRA(1) = [24.2737367272834];
% X0 = [0.0166290619855227;0.693438634162518;31.6912997059983];   %10Ail_2
% TRA(2) = [26.0435471267202];
% X0 = [0.0205060715159807;-0.202673986703291;36.1988955357060];  %10Ail_3
% TRA(3) = [21.9246736973074];
% X0 = [-0.0467565067470815;1.73649459920815;27.7849809892100];   %10Ail_4
% TRA(4) = [22.1307199973908];
% X0 = [-0.174562455198627;5.73003930980793;5.40242627289988];    %10Ail_5
% TRA(5) = [21.4454024842013];

% X0 = [-0.176667630790858;5.79757784679388;5.01794033444710];
% % TRA(1)= [21.4395174980359]
% X0 = [-0.000974739431267007;0.798719459083451;27.2929166532464];
% TRA(1)= [21.6058160404964]
% X0 = [-0.0115660431162276;0.922989896393046;28.4125750597272];
% TRA(1)= [21.5390873880433]
% X0 = [-0.0200321024205869;1.80499398865993;25.3629778699368];
% TRA(1)= [25.6229123042694];
% X0 = [-0.176667630790858;5.79757784679388;5.01794033444710];
% TRA(1)= [24.4304217974120];


%% SIMULATE GUST RESPONSE OF ORIGINAL DESIGN
main_pre_MF03_AilVar
main_processing
[u,x,y] = main_LQR(ssLinRed,StateWORBRed,SimDef,GustDef);
zDispWORBRed = calc_disprot(u,x,y,BeamDef,SimDef,trimstate,StateWORBRed,StateWORB,Nmodes,BeamSeg,Mode,Ti,RVnum);

% Calculate total wing mass
BaseRA = max(abs(BeamProp.RA));
% Calculate strain energy of initial configuration
maxMomOri = max(abs(y(1,:)));


%% OPTIMISATION
[history,searchdir] = runfmincon_MF03_TM_MMax_AilVar(maxMomOri,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,X0,BaseRA);
timeSimulation = toc;
filename = 'MF03_TM_10MMax_10Ail_2_new.mat';
save(filename);
