%% DEFINE DESIGN VECTOR
% Design vector contains three coefficient[a;b;c] of quadratic mass
% distribution  RA(i) = ax(i)^2+bx(i)+c
tic
X0 = [0;0;35.71] ;
% X0 = [-1.444428e-02;4.401484e-01;3.513538e+01]; % 5% MMax Reduction
% X0 = [-1.562341e-01;5.041303e+00;8.271621e+00]; % 10% MMax Reduction
% X0 = [-0.0627004898;1.8228122943;35.592814046]; % 20% MMax Reduction
% X0 = [-1.084040e-01;3.351945e+00;3.190871e+01]; % 30% MMax Reduction
% X0 = [-1.612312e-01;4.624077e+00;3.689332e+01]; % 40% MMax Reduction
% X0 = [-3.100421e-01;9.117161e+00;2.580147e+01]; % 50% MMax Reduction


% X0 = [-0.00484600708860643;0.442227438085148;34.1365244670021];   % 10% Span aileron
% X0 = [-0.00935383739247140;0.455372158884398;35.4794231714968];   % 20% Span aileron
% X0 = [-0.0149687225168727;0.568416310470364;35.0092321022422]     % 30% Span aileron
% X0 = [-0.177883325760070;5.74002098448955;5.00001780771681]       % 40% Span aileron
% X0 = [-0.175855942272830;5.69206453315224;5.00076443210118]       % 50% Span aileron
%% SIMULATE GUST RESPONSE OF ORIGINAL DESIGN
main_pre_MF03
main_processing
[u,x,y] = main_LQR(ssLinRed,StateWORBRed,SimDef,GustDef);
zDispWORBRed = calc_disprot(u,x,y,BeamDef,SimDef,trimstate,StateWORBRed,StateWORB,Nmodes,BeamSeg,Mode,Ti,RVnum);

% Calculate total wing mass
BaseRA = max(abs(BeamProp.RA));
% Calculate strain energy of initial configuration
maxMomOri = max(abs(y(1,:)));


%% OPTIMISATION
[history,searchdir] = runfmincon_MF03_TM_MMax(maxMomOri,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,X0,BaseRA);
timeSimulation = toc;
filename = 'MF03_TM_10MMax_20Ail.mat';
save(filename);
