function [StructMat,massmats1] = beamconvergence(BeamDef,BeamProp)
% This function calculates the natural displacement modes of the wing. 
% Assume there is no external force and linearise about the unloaded 
% equilibrium of x1 = x2 = 0. Refer to Eqn. (3.1). 
% This natural modes will be used to describe the aeroservoelastic
% system in modal coordinates(Eqn. 3.55).
%
% Solve : 
% mx1 - x2' - ex2  = 0
% cx2 - x1' + etx1 = 0
%
% Using separation variable (x1 = sin, x2 = cos)
% 
% ==> wm phi1 = d(phi2)/ds + E phi2
%    -wc phi2 = d(phi1)/ds - E^T phi1
% 
% (mw) and (-cw) in one matrix. d/ds + E ; d/ds -E^T in the other.
% 
% The result is arranged as x1;x2;x1;x2
% Child coordinate system as column vectors in the frame of the parent
%
% Input : 
    % BeamDef  : Beam Definition
        % NumNode   : Number of nodes
        % Ndivs     : Number of division
        % Ndl       : Length of element
        % Length    : Length of beam
    % BeamProp  : Beam properties  
        % RA        : Mass per unit length
        % RI1       : Sectional moment inertia per unit length(x direction)
        % RI2       : Sectional moment inertia per unit length(y direction)
        % RI3       : Sectional moment inertia per unit length(z direction)
        % EA        : Young Modulus in x,y,z direction
        % EI1       : In-plane bending stiffness
        % EI2       : Out-plane bending stiffness 
        % GJ        : Twisting stiffness
        

% Output :
    % StructMat : Structural Matrix
        % Dreal     : Eigen value of structure(Natural frequency)
        % Vreal     : Eigen vector (Mode shape)
        % Kmat      : K matrix of wing
        % Mmat      : M matrix of wing

NumNode = BeamDef.NumNode;
Ndivs   = BeamDef.Ndivs;
Length  = BeamDef.Length;
Ndl     = BeamDef.Ndl;
RA      = BeamProp.RA;
RI1     = BeamProp.RI1;
RI2     = BeamProp.RI2;
RI3     = BeamProp.RI3;
EA      = BeamProp.EA;
EI1     = BeamProp.EI1;
EI2     = BeamProp.EI2;
GJ      = BeamProp.GJ;

StandardM = zeros(6,6,Ndivs);
StandardC = zeros(6,6,Ndivs);
for i = 1:Ndivs
    StandardM(:,:,i) = diag([RA(i),RA(i),RA(i),RI1(i),RI2(i),RI3(i)]);       % Mass matrix per unit length
    StandardC(:,:,i) = inv(diag([EA(i),EA(i),EA(i),GJ(i),EI2(i),EI1(i)]));   % Compliance matrix per unit length
    
    % For Goland wing setup. 
    % Since the c.g. is not the same with e.a. position, there is a coupling between twisting and bending
    StandardM(3,4,i) = 0.2*0.9144*RA(i);
    StandardM(4,3,i) = 0.2*0.9144*RA(i);
end


% Put basic properties of beam into 'seg' structure
seg(1).Nseg  = Ndivs;               % Number of division
seg(1).L     = Length;              % Length of wing
seg(1).M     = StandardM;           % Mass matrix per unit length
seg(1).C     = StandardC;           % Compliance matrix per unit length    
seg(1).Ncon  = 0;                   % Number of segments this segment connects to (with the end node, only)
seg(1).ConID = [];                  % List of ID of the child segment
Nbeam = 1;

% Reclaim original shape            
seg(1).ownC = eye(3);
for segID = 1:Nbeam
    for i = 1:seg(segID).Ncon
        seg(seg(segID).ConID(i)).ownC = seg(segID).ownC*seg(segID).ConMat(:,:,i);
    end
end

seg(1).sID      = 1;
seg(1).startID  = 1;
seg(1).dl       = seg(1).L/seg(1).Nseg;

curID       = 1;
curID       = curID+1;

for i = 1:Nbeam
    for j = 2:seg(i).Nseg+1
        seg(i).sID = [seg(i).sID;curID];
        curID = curID+1;
    end
    % Fill the start of all connected segments          %?????????????%
    for j = 1:seg(i).Ncon
        seg(seg(i).ConID(j)).sID = curID-1;
    end
end


%============================= Define Kmat ===============================% 
% Nodal force = Kmat * displacement
% Use 'beamsectioncoeff' to extract this information

Kmat  = zeros(NumNode*6,NumNode*6);         % Kmat initialisation
segID = 1;
    
% First obtain transformation matrix
% ownC*local = global
segC = Rot6(seg(segID).ownC);
% each segment is an element
% Kmat is ordered by node, *6
for i = 1:seg(segID).Nseg
    % for i=1:1
    % start and end of beam element
    baseidstart = (seg(segID).sID(i)-1)*6;
    baseidend   = (seg(segID).sID(i+1)-1)*6;
    [fstartmat,fendmat,fstartmatr,fendmatr] = beamsectioncoeff(seg(segID).dl,seg(segID).C(:,:,i));
    % start and end converts a set displacement at the end of the
    % segment to the equivalent forces

    % first, displacement at the end node
    % localforce  = fendmat*localdisplacement
    % globalforce = ownc*localforce
    % localdisp   = ownc'*globaldisp
    fendmatg    = segC*fendmat*segC';
    fstartmatg  = segC*fstartmat*segC';
    Kmat(baseidend+1:baseidend+6,baseidend+1:baseidend+6)     = Kmat(baseidend+1:baseidend+6,baseidend+1:baseidend+6)+fendmatg;
    Kmat(baseidstart+1:baseidstart+6,baseidend+1:baseidend+6) = Kmat(baseidstart+1:baseidstart+6,baseidend+1:baseidend+6)+fstartmatg;

    % then, displacement at the start node
    % Follows same process
    fendmatgr   = segC*fendmatr*segC';
    fstartmatgr = segC*fstartmatr*segC';
    Kmat(baseidend+1:baseidend+6,baseidstart+1:baseidstart+6)     = Kmat(baseidend+1:baseidend+6,baseidstart+1:baseidstart+6)+fendmatgr;
    Kmat(baseidstart+1:baseidstart+6,baseidstart+1:baseidstart+6) = Kmat(baseidstart+1:baseidstart+6,baseidstart+1:baseidstart+6)+fstartmatgr;
end
%=========================================================================%   

%============================= Define Mmat ===============================%
% Initialisation of Mass Matrix
NDof = NumNode*6;               
Mmat = sparse(NDof,NDof);

% Calculte total mass matrix per wing element
massmats1(:,:,1:Ndivs) = StandardM(:,:,1:Ndivs)*Ndl;

% Direct Mass Lumping
% Mass matrix at the first node (root)
Mmat(1:6,1:6) = massmats1(:,:,1)/2;

for i = 2:Ndivs
    Mmat((i-1)*6+1:i*6,(i-1)*6+1:i*6) = massmats1(:,:,i-1)/2+massmats1(:,:,i)/2;
end

% Mass matrix at the last node (tip)
i = Ndivs+1;
Mmat((i-1)*6+1:i*6,(i-1)*6+1:i*6) = massmats1(:,:,Ndivs)/2;

% Assembly mass matrix
Mmat = full(Mmat);
%=========================================================================%   



%========================== Calculate eigenvector ========================%   
% w2Mmat=Kmat
[EVecs,EVals] = eig(Kmat(7:end,7:end),Mmat(7:end,7:end));       
EVecs = [zeros(6,size(Mmat,1)-6);EVecs];

% Compare in-plane part of the modes with out-of-plane part of the modes
keepmode = [];
for i = 1:size(EVals,1)
    EV = reshape(EVecs(:,i),6,Ndivs+1);
    EVs = sum(abs(EV),2);
    winplane = EVs(2) + EVs(6);
    wothers = EVs(1) + EVs(3) + EVs(4) + EVs(5);
    if winplane < wothers
        keepmode = [keepmode,i];
    end
end

EVecs = EVecs(:,keepmode);
EVals = EVals(keepmode,keepmode);

Dreal = diag(EVals);
Vreal = EVecs;

StructMat.Dreal = Dreal;
StructMat.Vreal = Vreal;
StructMat.Kmat = Kmat;
StructMat.Mmat = Mmat;

return              
