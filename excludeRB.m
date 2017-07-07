function [wingSSG,StateWORB] = excludeRB(MC,EigMatFull,SBG,SCG,trimstate,RVnum,Nmodes)

AMat = MC.AMat;

RV = 1:RVnum;   % Vector index exclude rigid body state
ReductionVec  = [RV,RV+Nmodes+6,RV+2*Nmodes+6,RV+2*Nmodes+6+Nmodes+6];  % Pick RVnum for each Velocity,Moment,and Aerodynamics state
ReductionVecS = [RV,RV+Nmodes+6];
AMatWORB = AMat(ReductionVecS,ReductionVecS);        % Matrix A without element related to rigid body state
EigMatWORB = EigMatFull(ReductionVec,ReductionVec);  % Pre-dynamic matrix without element related to rigid body states
SBGWORB = SBG(ReductionVec,:);                       % Input matrix according to reduction vector
SCGWORB = SCG(:,ReductionVec);                       % Output matrix according to reduction vector
StateWORB = trimstate(ReductionVec);                 % State 

EigMatWORB(1:numel(ReductionVecS),:) = AMatWORB\EigMatWORB(1:numel(ReductionVecS),:);   % Pre-dynamic matrix multiply by AMatRed
SBGWORB(1:numel(ReductionVecS),:) = AMatWORB\SBGWORB(1:numel(ReductionVecS),:);             % Sub Block of Input matrix multiply by AMatRed

%% STATE SPACE 
wingSSG = ss(EigMatWORB,SBGWORB,SCGWORB,zeros(size(SCGWORB,1),size(SBGWORB,2)));

end
