 function [InputMatrix,OutputMatrix] = systemdxioGust(MC,BeamDef,GustDef,AeDef,BeamSeg,Mode,trimstate,Nmodes)
%% THIS FUNCTION CALCULATES INPUT AND OUTPUT MATRIX FOR GUST RESPONSE
% INPUT :
    % Ly            & 1x1   & Gust spatial characteristics 
    % trimstate     & 2*Nmodes+6+nae*(Nmodes+6) x 1  & Trim State 
    % NumNode       & 1x1   & Number of Node
    % neig          & 1x1   & Number of retained mode
    % nae           & 1x1   & Number of aerodynamic states
    % Nflaps        % 1x1   & Number of flap
    % Mode          & struct& Mode database
    % moneord       & Nmodes+6 x Nmodes+6 x Nmodes+6 & Matrix H1
    % mtwoord       & Nmodes+6 x Nmodes+6 x Nmodes+6 & Matrix H2
    % mthreeord     & Nmodes+6 x Nmodes+6 x Nmodes+6 x Nflaps & Matrix H3
    % trimFlap      & 1x1   & Trim deflection of flap
    % Vi            & 1x1   & Freestream velocity
    % aeAb          & 1xnae & A multiply by b(theodorsen solution)
    % BeamSeg       & struct& Beam properties database

% OUTPUTS:
    % InputMatrix  & 2*Nmodes+6+nae*(Nmodes+6) x (Nflaps+2) & Input matrix
    % OutputMatrix & 3 x 2*Nmodes+6+nae*(Nmodes+6) & Output matrix
    
NumNode = BeamDef.NumNode;
Ly      = GustDef.Ly;
nae     = AeDef.nae;
Vi      = AeDef.Vi;
aeAb    = AeDef.aeAb;
Nflaps  = AeDef.Nflaps;

moneord     = MC.moneord;
mtwoord     = MC.mtwoord;
mthreeord   = MC.mthreeord;

x = trimstate(1:2*Nmodes+6);

% Velocity-Angular velocity state
x1 = x(1:Nmodes+6);

% Aerodynamic state
xae = trimstate(2*Nmodes+6+1:2*Nmodes+6+nae*(Nmodes+6));
% Convert ae states to matrix form for better manipulation
xam = zeros(Nmodes+6,nae*(Nmodes+6));
for i = 1:nae
    xam(:,(Nmodes+6)*(i-1)+1:(Nmodes+6)*i) = eye(Nmodes+6)*aeAb(i);
end
xaa = xam*xae;

%% COMPUTE INPUT MATRIX
% FLAP
dFlap = zeros(2*Nmodes+6+nae*(Nmodes+6),Nflaps);

for k = 1:Nmodes+6
    for nf = 1:Nflaps
        dFlap(k,nf) = dFlap(k,nf) + x1'*mthreeord(:,:,k,nf)*x1;
    end
end

% THRUST
% Compute unity thrust vector
UnitThrust = [0;1;0;0;0;0];
synthrust = zeros(2*Nmodes+6+nae*(Nmodes+6),1);
for i = 1:NumNode
    for j = 1:Nmodes+6
        synthrust(j) = synthrust(j)+dot(UnitThrust,Mode(j).Phi1(:,i));
    end
end

% GUST EFFECT 
% Spatial distribution
SpDist = zeros(3,numel(NumNode));
for i = 1:NumNode
%     SpDist(3,i) = 1/2*cos(pi*BeamSeg.NodeL(i)/Ly);    % Non uniform spatial distribution 
    SpDist(3,i) = 1;        % Uniform spatial distribution
end

Q1g = zeros(Nmodes+6,1);
for j = 1:Nmodes+6
    for s = 1:NumNode
        for i = 1:NumNode
            Q1g(j,1) =  Q1g(j,1) + Mode(j).Phi1(1:3,s)'*SpDist(:,i);
        end
    end
end

% Due to aerodynamics 
SB_q1 =  zeros(Nmodes+6,Nmodes+6);
for j = 1:Nmodes+6
        SB_q1(j,1:Nmodes+6) = x1'*moneord(:,:,j)+x1'*moneord(:,:,j)'+(mtwoord(:,:,j)*(xaa*Vi))';
end

dGust = zeros(2*Nmodes+6+nae*(Nmodes+6),1);
dGust(1:numel(Q1g)) = SB_q1*(-Q1g);

%% SINCE WE HAVE THE TRIM POSITION OF FLAP IS ZERO (TRIMFLAP = 0), WE DO NOT NEED TO CALCULATE THIS TERMS
% for k = 1:Nmodes+6
%     dGust(1:Nmodes+6,1) = dGust(1:Nmodes+6,1) + (x1'*moneord(:,:,k))' + (x1'*moneord(:,:,k)')' + (mtwoord(:,:,k)*(xaa*Vi));
% end                       

% Due to flap deflection
% for k = 1:neig+6
%     for nf = 1:Nflaps
%         dGust(1:neig+6,1) = dGust(1:neig+6,1) + (x1'*mthreeord(:,:,k,nf))'*trimFlap + (x1'*mthreeord(:,:,k,nf)')'*trimFlap;
%     end
% end

%% ASSEMBLY INPUT MATRIX
InputMatrix = [dFlap,synthrust,dGust];

%% OUTPUT MATRIX
OutputMatrix = zeros(3,2*Nmodes+6+nae*(Nmodes+6));      % Initialisation

% Populate output matrix
for i = 1:Nmodes
    % Moment
    OutputMatrix(1,i+Nmodes+6) = OutputMatrix(1,i+Nmodes+6)+Mode(i).Seg(1).phi2(5,1);
end

for i = 1:Nmodes+6
    % Velocity
    OutputMatrix(2,i) = OutputMatrix(1,i)+Mode(i).Phi1(3,end);
end

for i = 1:Nmodes
    % Force
    OutputMatrix(3,i+Nmodes+6) = OutputMatrix(3,i+Nmodes+6)+Mode(i).Seg(1).phi2(3,1);
end
