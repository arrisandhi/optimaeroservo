function  [ModalRE, ModalCE, dRH, dCH, dR, dC] = dispintegrationvar_Arri(CurrentState, nInterest, ttime, NumNode, neig, NumIndices, BeamSeg, Nseg, Mode)
% This function returns diplacement along the wing according to the
% CurrentState and obtain the variation of displacement with the next state
% INPUT :
%   CurrentState  : Current State
%   NumNode       : Number of node
%   neig          : Number of modes to be retained
%   NumIndices    : Number of actual indices, including duplication
%   BeamSeg       : Beam segment data structure
%   nInterest     : Node interest
%   ttime         : Time to be considered
%
% OUTPUT :
%   ModalRE   : Diplacement along the wing according to CurrentState
%   ModalCE   : Curvature along the wing according to CurrentState
%   dRH       : Displacement variation
%   dCH       : Curvature variation

Y = CurrentState';

% cycle through each strain mode
dRH = zeros(3,NumNode,neig);
dCH = zeros(3,3,NumNode,neig);
for k = 1:neig
    dR = zeros(3,NumNode);
    dC = zeros(3,3,NumNode);
    ModalRE = zeros(3,NumNode);
    ModalCE = zeros(3,3,NumNode);
    Yvar = zeros(size(Y));
    Yvar(k+neig+6) = 1;
    Vs = zeros(6,NumIndices);
    dVs = zeros(6,NumIndices);
    
    % Obtain collocated displacement derivative vector. Refer to Eqn. 3.47 in dissertation.
    for i = 1:neig
        for segi = 1:Nseg
            Vs(:,BeamSeg(segi).NodeStart:BeamSeg(segi).NodeStart+BeamSeg(segi).NodeNo-2) = Vs(:,BeamSeg(segi).NodeStart:BeamSeg(segi).NodeStart+BeamSeg(segi).NodeNo-2)+Mode(i).Seg(segi).cphi2*Y(i+neig+6);
            dVs(:,BeamSeg(segi).NodeStart:BeamSeg(segi).NodeStart+BeamSeg(segi).NodeNo-2) = dVs(:,BeamSeg(segi).NodeStart:BeamSeg(segi).NodeStart+BeamSeg(segi).NodeNo-2)+Mode(i).Seg(segi).cphi2*Yvar(i+neig+6);
        end
    end
    
    
    for segi = 1:Nseg
        % start from first node of each segment
        j = 1;
        % Except starting segment, whose rotation is already defined
        if segi > 1
            % Transform coordinates into the current beam's starting coords
            % Defined still in global frame
            StartCoord = ModalCE(:,:,BeamSeg(BeamSeg(segi).Parent).NodeOrder(end));
            StartR = ModalRE(:,BeamSeg(BeamSeg(segi).Parent).NodeOrder(end));
            ModalRE(:,BeamSeg(segi).NodeOrder(1),ttime) = StartR;
            ModalCE(:,:,BeamSeg(segi).NodeOrder(1),ttime) = StartCoord; 


            dR0 = dR(:,BeamSeg(BeamSeg(segi).Parent).NodeOrder(end));
            dC0 = dC(:,:,BeamSeg(BeamSeg(segi).Parent).NodeOrder(end));
            dR(:,BeamSeg(segi).NodeOrder(1)) = dR0;
            dC(:,:,BeamSeg(segi).NodeOrder(1)) = dC0;

        else
            ModalRE(:,1,ttime) = [0;0;0];
            ModalCE(:,:,1,ttime) = eye(3);
            StartR = [0;0;0];
            StartCoord = eye(3);

            dR0 = [0;0;0];
            dC0 = zeros(3,3);

            dR(:,1) = dR0;
            dC(:,:,1) = dC0;

        end

        % Find running direction, transform 1,0,0 into GLOBAL AXES
        edir = BeamSeg(segi).GlobalAxes*[1;0;0];

        % Now integrate along beam
        for j = 2:BeamSeg(segi).NodeNo
            % Extract g0 and k0           
            % REMEMBER GAMMA IS COLLOCATED, THUS DIVIDE BY DL
            g0 = Vs(1:3,BeamSeg(segi).NodeStart-2+j)/BeamSeg(segi).NodeDL(j-1);
            k0 = Vs(4:6,BeamSeg(segi).NodeStart-2+j)/BeamSeg(segi).NodeDL(j-1);
            % Transform strains and curvatures too.
            % Different from angular velocity, as strains are defined in local
            % coords at the beginning
            % But as the integration occurs in global
            tg0 = BeamSeg(segi).GlobalAxes*g0;
            tk0 = BeamSeg(segi).GlobalAxes*k0;

            dg0 = dVs(1:3,BeamSeg(segi).NodeStart-2+j)/BeamSeg(segi).NodeDL(j-1);
            dk0 = dVs(4:6,BeamSeg(segi).NodeStart-2+j)/BeamSeg(segi).NodeDL(j-1);
            dgamma = BeamSeg(segi).GlobalAxes*dg0;
            dkappa = BeamSeg(segi).GlobalAxes*dk0;

            [StartR,StartCoord,dR0,dC0] = rkelesvar(StartR,StartCoord,tg0,tk0,edir,BeamSeg(segi).NodeDL(j-1),dR0,dC0,dgamma,dkappa);
            ModalRE(:,BeamSeg(segi).NodeOrder(j),ttime) = StartR;
            ModalCE(:,:,BeamSeg(segi).NodeOrder(j),ttime) = StartCoord;

            dR(:,BeamSeg(segi).NodeOrder(j)) = dR0;
            dC(:,:,BeamSeg(segi).NodeOrder(j)) = dC0;
        end
    end
% pull-back operation into reference node
    for j=1:NumNode
        dCH(:,:,j,k) = ModalCE(:,:,nInterest,ttime)'*dC(:,:,j)+dC(:,:,nInterest)'*ModalCE(:,:,j,ttime);
        dRH(:,j,k) = ModalCE(:,:,nInterest,ttime)'*(dR(:,j)-dR(:,nInterest))+dC(:,:,nInterest)'*(ModalRE(:,j,ttime)-ModalRE(:,nInterest,ttime));
    end
end


MCE = ModalCE(:,:,nInterest);
for i = 1:NumNode
    ModalCE(:,:,i) = MCE'*ModalCE(:,:,i);
end

% Now dCH and dRH are variations from trim condition (assumed to be Y) as a
% function of x2



% Gr=C(:,:,nInterest)'*gvector

% <(ModalCE+dCH*x2)'*Gr,MPhi1>
% = MPhi1'*(ModalCE+dCH*x2)'*Gr

global Grav0
global Grav1
Grav0 = zeros(neig+6,3);
Grav1 = zeros(neig+6,neig,3);
for i = 1:neig+6
    for j = 1:NumNode
        Grav0(i,:) = Grav0(i,:)+Mode(i).MPhi1(1:3,j)'*ModalCE(:,:,j)';
        for k = 1:neig
            Grav1(i,k,:) = Grav1(i,k,:)+reshape(Mode(i).MPhi1(1:3,j)'*dCH(:,:,j,k)',1,1,3);
        end
    end
end

