%% FOR THE PURPOSE OF MODEL ORDER REDUCTION(MOR) PROCESS, SEPARATE EIGMATLIN AND EIGMATNONLIN
EigMatLinWORB = EigMatLin(ReductionVec,ReductionVec);
EigMatNonLinWORB = EigMatNonLin(ReductionVec,ReductionVec);
EigMatLinWORB(1:numel(ReductionVecS),:) = AMatWORB\EigMatLinWORB(1:numel(ReductionVecS),:);
EigMatNonLinWORB(1:numel(ReductionVecS),:) = AMatWORB\EigMatNonLinWORB(1:numel(ReductionVecS),:);

% Linear part of state space
wingSSGLinMOR = ss(EigMatLinWORB,SBGWORB,SCGWORB,zeros(size(SCGWORB,1),size(SBGWORB,2)));

%% 
% EigMatLin(1:2*neig+6,1:2*neig+6) = AMat\EigMatLin(1:2*neig+6,1:2*neig+6);
% 
% wingSSGLinMORFull = ss(EigMatLin,SBG,SCG,zeros(size(SCG,1),size(SBG,2)));

%% CREATE BALANCE STATE SPACE FROM LINEAR PART OF STATE SPACE
[ssLinBal,g,Tr,Ti] = balreal(wingSSGLinMOR);
k = find(g<1e-6)
TrMod = zeros(size(Tr));
TiMod = TrMod;
TrMod(1:k(1)-1,1:k(1)-1) = Tr(1:k(1)-1,1:k(1)-1);
TiMod(1:k(1)-1,1:k(1)-1) = Ti(1:k(1)-1,1:k(1)-1);
%% MODEL REDUCTION
% CALCULATE TRANSFORMED NONLINEAR TERMS
Qr = TrMod*EigMatNonLinWORB*TiMod*TiMod;
Qr = Qr(1:k(1)-1,1:k(1)-1);

% Determine state reduction parameter
ssLinRed = modred(ssLinBal,k,'Truncate');

%% % RRi Tr A Ti RR, RRi Tr, Ti RR
% [ssLinRed,v1,s1,t1] = slowO(ssLinRed,k(1)-1);

% redsize = size(ssLinRed1.a,1);
% Apre = v1'*RRi*Tr;
% Apost = Ti*RR*v1;

% Determine state reduction parameter
% k = find(g<1e-6);
% redorder = 120;
% ssLinRed2 = modred(ssLinBal,redorder+1:msize,'Truncate');
% 
% % RRi Tr A Ti RR, RRi Tr, Ti RR
% [ssLinRed2,v12,s12,t12] = slowO(ssLinRed2,k(1)-1);
% 
% redsize = size(ssLinRed2.a,1);
% Apre = v12'*RRi*Tr;
% Apost = Ti*RR*v12;

%%