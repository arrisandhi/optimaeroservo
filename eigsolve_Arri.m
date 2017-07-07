function [goneord, gtwoord, moneord, mtwoord, mthreeord, dcoefford, EigMat, wingSS, SB, SC, trimstate] = eigsolve_Arri(Vi, neig, nae, btoV, aeomega, gamma1, gamma2, muone, mutwo, muthree, Nflaps, AMat, BMat, aeAb, NumNode, Mode, RVnum)
%% This function returns trimstate, dynamic matrix, input matrix, and output
% matrix of aeroservoleastic system.
% INPUT :
% Vi	 	& 1x1                            & Freestream velocity  
% neig		& 1x1                            & Number of modes 
% nae		& 1x1                            & Number of aerodynamic state
% btoV		& 6x6                            & Transformation matrix from rb to velocity 
% aeomega	& nae x 1                        & aebvector divided by aebvalue 
% alpha1 	& Nmodes+6 x Nmodes+6            & Matrix A1 
% alpha2  	& Nmodes x Nmodes                & Matrix A2 
% beta1     & Nmodes+6 x Nmodes              & Matrix WD
% beta2     & Nmodes                         & Nmodes+6 -WD 
% gamma1	& Nmodes+6 x Nmodes+6 x Nmodes+6 & Matrix gamma1 
% gamma2	& Nmodes+6 x Nmodes x Nmodes     & Matrix gamma2 matrix 
% muone     & Nmodes+6 x Nmodes+6 x Nmodes+6 & Matrix H1-flap 
% mutwo   	& Nmodes+6 x Nmodes+6 x Nmodes+6 & Matrix H2-flap  
% muthree	& Nmodes+6 x Nmodes+6 x Nmodes+6 x Nflaps & Matrix H3-flap 
% Nflaps	& 1x1                            & Number of flap segment  
% AMat	 	& Nmodes+6+Nmodes x Nmodes+6+Nmodes     & Full Matrix A 
% BMat	 	& Nmodes+6+Nmodes x Nmodes+6+Nmodes     & Full Matrix WD 
% NumNode 	& 1x1           & Number of node 
% aeAb	 	& 2x1           & A multiply by b (Theodorsen solution)
% Mode		& 1x26 struct   & Mode shape 
% RVNum     & 1x1           & Number of modes retained 
%
% OUTPUT :
% goneord	& Nmodes+6 x Nmodes+6 x Nmodes+6    & Matrix gamma1 
% gtwoord	& Nmodes x Nmodes x Nmodes+6        & Matrix gamma2 matrix 
% moneord	& Nmodes+6 x Nmodes+6 x Nmodes+6    & Matrix H1-flap 
% mtwoord 	& Nmodes+6 x Nmodes+6 x Nmodes+6    & Matrix H2-flap  
% mthreeord	& Nmodes+6 x Nmodes+6 x Nmodes+6 x Nflaps & Matrix H3-flap 	
% dcoefford &  & 
% EigMat	& Nmodes x Nmodes    & Dynamic matrix 
% SB		& Nmodes x 3         & Input matrix
% SC 		& 2xnumel(trimstate) & Vector input 

%% Assembly trimstate
% trimstate        1 to neig    : velocity(velocity and angular velocity) state structural
% trimstate   neig+1 to neig+6  : velocity(velocity and angular velocity) state rigid body
% trimstate   neig+7 to 2*neig+6: force moment state
% trimstate 2*neig+7 to  end    : aerodynamics state

% associated with each velocity mode
trimstate = zeros(2*neig+6+nae*(neig+6),1);

% trimsate(neig+2,1) : velocity mode rigid body in y direction 
% btoV(2,2) : velocity mode in y direction
trimstate(neig+2,1) = Vi/btoV(2,2);

% Definition of velocity state
x1 = trimstate(1:neig+6);

% Initialisation of aerodynamic state(written in nae column)
xaemat = zeros(neig+6,nae);

% Write Back AE STATES
for k = 1:neig+6
    for i = 1:nae
        xaemat(k,i) = x1(k)/Vi/aeomega(i,1);
    end
end

% Reshape aerodynamic state matrix and put them into trimstate vector
xae = reshape(xaemat,nae*(neig+6),1);
trimstate(2*neig+6+1:2*neig+6+nae*(neig+6)) = xae;




%% Assembly EigMat(NOT INCLUDING WD MATRIX AND HAS NOT BEEN MULTIPLIED BY AMAT_INVERSE)
EigMat = systemdxeigconvergence_Arri( trimstate, neig, nae, aeomega, aeAb, goneord, moneord, mtwoord, dcoefford, btoV);

%% Assembly EigMat. WD MATRIX TRANSFORMED INTO LAMBDAMAT. THEN ADDED INTO EIGMAT. THIS EIGMAT STILL HAS NOT BEEN CALCULATED BY AMAT_INVERSE
LambdaMat = zeros(2*neig+6+nae*(neig+6),2*neig+6+nae*(neig+6));
LambdaMat(1:2*neig+6,1:2*neig+6) = AMat*BMat; % Obtain Lambda matrix
EigMat = LambdaMat + EigMat;

%% Assembly EigMat. Effect of gust


%% For 
% Modes to retain
RV = 1:RVnum;
ReductionVec  = [RV,RV+neig+6,RV+2*neig+6,RV+2*neig+6+neig+6];
ReductionVecS = [RV,RV+neig+6];
AMatRed = AMat(ReductionVecS,ReductionVecS);
EigMat = EigMat(ReductionVec,ReductionVec);

% First block of EigMat needs multiplying by A inverse
EigMat(1:numel(ReductionVecS),:) = AMatRed\EigMat(1:numel(ReductionVecS),:);


% Obtain input and output matrix
[SB,SC] = systemdxioall_Arri( trimstate, neig, nae, mthreeord, NumNode, Mode, Nflaps);
SB = SB(ReductionVec,:);
SC = SC(:,ReductionVec);

% Input matrix needs multiplying by A inverse
SB(1:numel(ReductionVecS),:) = AMatRed\SB(1:numel(ReductionVecS),:);

% Cutoff frequency in rad/s
fcut = 200;

wingSS = ss(EigMat,SB,SC,zeros(size(SC,1),size(SB,2)));
[wingSS,~] = freqsep(wingSS,fcut);
return