function [AMat,goneord,gthreeord,mthreeord] = guyanRecomp(BeamDef,BeamSeg,BeamProp,AeDef,Mode,Nmodes)
% Get ready for optimisation, recompute dynamic matrix of system based on a standard
% set of modes. To do this, we require to recompute alpha1/alpha2,
% gone/dcoeff, and mthree.
% INPUT :
% NumNode   & 1x1 & Number of node	
% Nmodes	& 1x1 & Number of modes
% neig		& 1x1 & Number of modes
% Mode		& 1x26struct 	& Modes 
% EAn		& 1x1 & Young modulus 
% EI1n		& 1x1 & Torsional modulus(z direction)
% EI2n		& 1x1 & Torsional modulus(y direction)
% GJn		& 1x1 & Torsional modulus(x direction)
% RAn		& 1x1 & Torsional modulus(x direction)
% RI1n		& 1x1 & Torsional modulus(x direction) 
% RI2n		& 1x1 & Torsional modulus(x direction) 
% RI3n		& 1x1 & Torsional modulus(x direction)
% BeamSeg   & struct& BeamSeg information
% Nflaps 	& 1x1 & Number of flap
% aeavalue  & 1x1 & Ratio between aerodynamic center and half chord	
% aebvalue  & 1x1 & Half chord length
% rhoinf	& 1x1 & Air density	
% CLD	    & 1xNflaps & Flap drag coefficient	
% CMD	    & 1xNflaps & Flap moment coefficient	
% flapdist  & Nflaps x NumNode & Flap distribution. Assign value 1 at node index for flap position
%
% OUTPUT:
% goneord	& Nmodes+6 x Nmodes+6 x Nmodes+6          & Matrix gamma1 \\ \hline
% dcoefford &                                         & Matrix gamma2
% mthreeord	& Nmodes+6 x Nmodes+6 x Nmodes+6 x Nflaps & Matrix H3-flap \\ \hline	
% AMat	 	& Nmodes+6+Nmodes x Nmodes+6+Nmodes       & Full Matrix A \\ \hline

NumNode = BeamDef.NumNode;
Ndivs = BeamDef.Ndivs;
Ndl = BeamDef.Ndl;
Nseg = BeamDef.Nseg;

EAn     = BeamProp.EA;
EI1n    = BeamProp.EI1;
EI2n    = BeamProp.EI2;
GJn     = BeamProp.GJ;
RAn     = BeamProp.RA;
RI1n    = BeamProp.RI1;
RI2n    = BeamProp.RI2;
RI3n    = BeamProp.RI3;

aeavalue    = AeDef.aeavalue;
aebvalue    = AeDef.aebvalue;
rhoinf      = AeDef.rhoinf;
CLD         = AeDef.CLD;
CMD         = AeDef.CMD;
Nflaps      = AeDef.Nflaps;
flapdistr   = AeDef.flapdistr;

StandardM = zeros(6,6,Ndivs);
StandardC = zeros(6,6,Ndivs);
for i = 1:Ndivs
    StandardM(:,:,i) = diag([RAn(i),RAn(i),RAn(i),RI1n(i),RI2n(i),RI3n(i)]);       % Mass matrix per unit length
    StandardC(:,:,i) = inv(diag([EAn(i),EAn(i),EAn(i),GJn(i),EI2n(i),EI1n(i)]));   % Compliance matrix per unit length %
end

% For Goland wing setup
% StandardM(3,4) = -0.2*0.9144*RAn;
% StandardM(4,3) = StandardM(3,4);

% Calculte total mass matrix per wing element
massmats1(:,:,1:Ndivs) = StandardM(:,:,1:Ndivs)*Ndl;

BeamM = zeros(6,6,NumNode);     % Lumped mass 
BeamC = zeros(6,6,NumNode);      % Lumped compliance here (C*dl)

BeamM(:,:,1) = massmats1(:,:,1)/2;
for N = 2:Ndivs
    BeamM(:,:,N) = massmats1(:,:,i-1)/2+massmats1(:,:,i)/2;
end
BeamM(:,:,NumNode) = massmats1(:,:,Ndivs)/2;

for N = 1:NumNode
    if N < NumNode
        BeamC(:,:,N) = StandardC(:,:,N)*BeamSeg(1).NodeDL(N);
    end
end

% StandardC = inv(diag([EAn(1),EAn(1),EAn(1),GJn(1),EI2n(1),EI1n(1)]));
% StandardM = diag([RAn(1),RAn(1),RAn(1),RI1n(1),RI2n(1),RI3n(1)]);
% % StandardM(3,4) = -0.2*0.9144*RAn;
% % StandardM(4,3) = StandardM(3,4);
% 
% for N = 1:NumNode
%     BeamM(:,:,N) = StandardM*BeamSeg(1).NodeDL(1);
%     if N == 1 || N == NumNode
%         BeamM(:,:,N) = BeamM(:,:,N)/2;
%     end
%     if N < NumNode
%         BeamC(:,:,N) = StandardC*BeamSeg(1).NodeDL(N);
%     end
% end


% Recompute alpha couplings
alpha1new = zeros(Nmodes+6,Nmodes+6);
alpha2new = zeros(Nmodes,Nmodes);
for j = 1:Nmodes+6
    % k=j;
    for k = 1:Nmodes+6
        alpha1new(j,k) = 0;
        if j <= Nmodes
            if k <= Nmodes
                alpha2new(j,k) = 0;                           
            end
        end
        % x1Mx1
        % In discrete form
        for N = 1:NumNode
            alpha1new(j,k) = alpha1new(j,k)+Mode(j).Phi1(:,N)'*Mode(k).MPhi1(:,N);          
        end
        for segi = 1:Nseg
            for N = 1:BeamSeg(segi).NodeNo-1           
                if j <= Nmodes && k <= Nmodes
                        % x2Cx2
                        alpha2new(j,k) = alpha2new(j,k)+Mode(j).Seg(segi).phi2(:,N)'*Mode(k).Seg(segi).cphi2(:,N);
                end
            end            
        end
        % [j k]
    end
end

AMat = blkdiag(alpha1new,alpha2new);



% Couplings
% Only goneord(:,neig+2,:) and goneord(neig+2,:,:) are needed since gamma1(k,i,j) == goneord(i,j,k)

gamma1 = zeros(Nmodes+6,Nmodes+6,Nmodes+6);
% Integrate x1L1x1
for N = 1:NumNode
    
    % gamma1
    gm1 = zeros(Nmodes+6,Nmodes+6,Nmodes+6);
    for k = 1:Nmodes+6
        for l = 1:Nmodes+6
            % pre-evaluate
            L1kx1l = Mode(k).L1P1(:,:,N) * Mode(l).MPhi1(:,N);
            for j = 1:Nmodes+6
                % x1L1(x1)Mx1
                % Discrete
                gm1(j,k,l) = Mode(j).Phi1(:,N)'*L1kx1l;
            end
        end
    end
    gamma1 = gamma1+gm1;
    % N
end



% Only dcoefford(neig+2,:,:) are needed because:
% dcoefford(i,j,k) == dcoeff(k,i,j);
% dcoeff(j,k,l) = gtwo(k,j,l)
% Thus, dcoefford(i,j,k) = gtwo(i,k,j))


% Global dcoefford
gthreeord = zeros(Nmodes+6,Nmodes,Nmodes);

% Gamma2 needed for deltacoeff
segi = 1;
for N = 1:BeamSeg(segi).NodeNo-1          
    dcf = zeros(Nmodes+6,Nmodes,Nmodes);
    for k = 1:Nmodes
        for l = 1:Nmodes    
            L2kx2l = Mode(k).L2P2(:,:,segi,N) * BeamC(:,:,N)*Mode(l).Seg(segi).phi2(:,N);
            j = Nmodes+2;
            
            % x1L2(x2)Cx2
            dcf(j,l,k) = Mode(j).Seg(segi).phi1(:,N)'*L2kx2l;
        end
    end
    gthreeord = gthreeord+dcf;
end

% Using flap definition to calculate matrix below
% mthreeord(i,j,k,nf) = mthree(k,i,j,nf);
% only mthreeord(:,neig+2,:) and mthreeord(neig+2,:,:) are needed

% Global mthreeord
mthreeord = zeros(Nmodes+6,Nmodes+6,Nmodes+6,Nflaps);

segi = 1;

% semi-DL
NodeDLh = ([BeamSeg(segi).NodeDL/2;0]+[0;BeamSeg(segi).NodeDL/2]);

CpPhi1 = zeros(6,Nmodes+6,BeamSeg(segi).NodeNo);
A3CpPhi1 = zeros(6,6,Nmodes+6,BeamSeg(segi).NodeNo,Nflaps);

% pre-compute transformed velocities and matrices
% This already exists
for k = 1:Nmodes+6
    for b = 1:BeamSeg(segi).NodeNo
        CpPhi1(:,k,b) = Mode(k).Phi1(:,BeamSeg(segi).NodeOrder(b));
        for nf = 1:Nflaps
            A3CpPhi1(:,:,k,b,nf) = A3funD(Mode(k).Phi1(:,BeamSeg(segi).NodeOrder(b)),CLD(nf),CMD(nf),aeavalue,aebvalue);
        end
    end
end

for i = 1:Nmodes+6
    for j = 1:Nmodes+6
        if (i==Nmodes+2) || (j==Nmodes+2)
            mdd = zeros(6,BeamSeg(segi).NodeNo,Nflaps);
            % Pre-compute second-order couplings
            for b = 1:BeamSeg(segi).NodeNo
                ndlhb = NodeDLh(b);
                % Multiple flaps
                for nf = 1:Nflaps
                    mdd(:,b,nf) = (A3CpPhi1(:,:,i,b,nf)*(Mode(j).Phi1(:,BeamSeg(segi).NodeOrder(b))))*ndlhb;
                end
            end
            for k = 1:Nmodes+6
                magd = zeros(Nflaps,1);

                for b = 1:BeamSeg(segi).NodeNo
                    cpp1p = CpPhi1(:,k,b)';

                    for nf = 1:Nflaps
                        magd(nf) = magd(nf)+cpp1p*mdd(:,b,nf)*flapdistr(nf,b);
                    end
                end

                for nf = 1:Nflaps
                    mthreeord(i,j,k,nf) = mthreeord(i,j,k,nf)+magd(nf);
                end
            end
        end
    end
end
mthreeord = mthreeord*rhoinf*aebvalue;

end
