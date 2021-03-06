function [cnlin, cnleq] = confun_MF03_TM_MMMS(X0,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,maxMomOri)
aRA = X0(1);
bRA = X0(2);
cRA = X0(3);
NodeL = BeamDef.NodeL;
Ndivs = BeamDef.Ndivs;
Ndl = BeamDef.Ndl;

RA = zeros(Ndivs,1);
for i = 1:Ndivs;
    RA(i) = aRA*NodeL(i)^2+bRA*NodeL(i)+cRA;
end

BeamProp.RA = RA;

%% CALCULATE STATIC MOMENT AS CONSTRAINT
wSeg = RA*Ndl*9.8;      %Calculate weight
MS = 0;
% Calculate moment
for i = 1:Ndivs
    MS = MS+ wSeg(i)*((i-1)*Ndl+Ndl/2);
end

if MS > 2e5
    cMS = 1;
    cMom = 1;
    cFS = 1;
else
    cMS = 0;
    %% OBTAIN K AND M MATRIX
    [StructMat] = beamconvergence(BeamDef,BeamProp);

    %% CALCULATE MATRIX COEFFICIENT
    [MC,~,AeDef,Mode,btoV] = guyancomp(Nmodes,BeamSeg,BeamDef,AeDef,StructMat);
    % [Vreal,alpha1,alpha2,beta1,beta2,goneord,gtwoord,gthreeord,moneord,mtwoord,mthreeord,aeomega,btoV,nae,neig,AMat,BMat,aeAb,Mode] = guyancomp(Nmodes, BeamSeg, Dreal, Vreal, Mmat, Nseg, NumIndices, NumNode, aeavalue, aebvalue, rhoinf, Nflaps, CLD, CMD,flapdistr);

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
    % [u,du,x,y] = main_MPCsptwir(ssLinRed,StateWORBRed,SimDef,GustDef);
    [u,x,y] = main_LQR(ssLinRed,StateWORBRed,SimDef,GustDef);
    maxMomRoot = max(abs(y(1,:)));
    
    % Nonlinear inequality constraints
    cMom = maxMomRoot - 0.80*maxMomOri;

    eigFlutter = max(real(eig(wingSSG.a)));
    if eigFlutter > 0
        cFS = 1;
    else
        cFS = 0;
    end
end

cnlin = [cMom ; cFS; cMS];
% Nonlinear equality constraints
cnleq = [];
end