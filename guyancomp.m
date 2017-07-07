function [MC,StructMat,AeDef,Mode,btoV] = guyancomp(Nmodes,BeamSeg,BeamDef,AeDef,StructMat)
%% This is the code to obtain the intrinsic model from a NASTRAN box beam, as in Palacios, Wang and Karpel, SDM2012.
% Input :
    % Nmodes        : The first number of modes that will be retained
    % BeamSeg       : Structure datatype contained properties of beam segment
    % Dreal         : Eigenvalue. Already in w^2 form
    % Vreal         : Eigenvector
    % Mmat          : Mass Matrix
    % Nseg          : Number of segment. Same with number of beam
    % NumIndices    : Number of actual indices, including duplications at connections
    % NumNode       : Number of nodes
    % aeavalue      : Ratio of distance of aerodynamic center(a.c) from elastic
    %                 axis(e.a) to half chord. If e.a. is backwards of a.c. then 
    %                 a is positive. The distance is a*b
    % aebvalue      : Half chord 
    % rhoinf        : Density of air
    % Nflaps        : Number of flap segment
    % CLD           : Drag coefficient of flap
    % CMD           : Moment coefficient of flap
    % flapdistr     : Flap position represented by value (1) at mass node

% Output :
    % Please see Eqn. 3.18, Eqn. 3.40 for detail definition
    % Vreal         : Eigenvector with size Nmodes+6
    % alpha1        : Matrix A1
    % alpha2        : Matrix A2
    % beta1         : Diagonal Matrix of Eigenvalue (W_D)
    % beta2         : Diagonal Matrix of Eigenvalue(minus sign) (W_D)
    % gamma1        : Matrix gamma1
    % gamma2        : Matrix gamma2
    % muone         : Matrix H1
    % mutwo         : Matrix H2
    % muthree       : Matrix H3
    % aeomega       : aebvector/aebvalue
    % btoV          : Transformation matrix from rb to velocity
    % nae           : Number of aerodynamics state
    % neig          : Number of considered eigenvalue.Equal to Number of modes
    % AMat          : Dynamics matrix of aeroelastic system
    % BMat          : Input matrix f aeroelastic system
    % aeAb          : Multiplication of aeAvector with aeBvector
    % Mode          : Mode shape. Contained variable below:
    %               Omega : Natural freq(From Dreal)
    %               Phi0  : Global displacement modes 
    %               Phi1  : Global velocity and angular velocity modes
    %               Phi2  : Global strain modes
    %               Phi0S : Split nodal Phi0
    %               Phi1S : Split nodal Phi1
    %               MPhi1 : Corresponding modes in momentum. Eqn. 5.19
    %                       Equal to integration of Psi1. Mass matrix times Phi1. 
    %                       Used for calculating A1 and gamma1 matrix
    %               For each segment contains:
    %                   cphi2 : Sectional curvature strains. Equal to integration of Psi2 
    %                   phi2  : local phi2
    %                   phi1  : local phi1, centerpoint averaged. 

Nseg        = BeamDef.Nseg;
NumIndices  = BeamDef.NumIndices;
NumNode     = BeamDef.NumNode;

aeavalue    = AeDef.aeavalue;
aebvalue    = AeDef.aebvalue;
rhoinf      = AeDef.rhoinf;
Nflaps      = AeDef.Nflaps;
CLD         = AeDef.CLD;
CMD         = AeDef.CMD;
flapdistr   = AeDef.flapdistr;

Dreal = StructMat.Dreal;
Vreal = StructMat.Vreal;
Mmat  = StructMat.Mmat;

% Dreal is the eigenvalues
% need w^2. This version's Dreal is already w^2

% Sort eigenvalue and eigenvector in ascending order
[Dreal,Ind] = sort(Dreal);
Vreal = Vreal(:,Ind);

% Take the first 'Nmodes' of modes
Dreal = Dreal(1:Nmodes);
Vreal = Vreal(:,1:Nmodes);

%=========================================================================%
%                     Rigid Body(RB) mode shapes                          %
%=========================================================================%
% This mode relates with forward flight motion
% Need one RB mode not originally in the system

% Make provisions for 6 rigid body modes. Put them at the last 6 columns
Nmodes = Nmodes+6;
Vreal  = [Vreal,zeros(size(Vreal,1),6)];
Dreal  = [Dreal;[0;0;0;0;0;0]];

%======================= REDEFINE RB MODES ===============================%

% Calculate centre of mass(xm)
xm = BeamSeg(1).NodeX(1,:)*full(Mmat(6*(BeamSeg(1).NodeOrder(1)-1)+1,6*(BeamSeg(1).NodeOrder(1)-1)+1));
sm = full(Mmat(6*(BeamSeg(1).NodeOrder(1)-1)+1,6*(BeamSeg(1).NodeOrder(1)-1)+1));
for segi = 1:Nseg
    for i = 2:BeamSeg(segi).NodeNo
        xm = xm + BeamSeg(segi).NodeX(i,:)*full(Mmat(6*(BeamSeg(segi).NodeOrder(i)-1)+1 , 6*(BeamSeg(segi).NodeOrder(i)-1)+1));
        sm = sm + full(Mmat(6*(BeamSeg(segi).NodeOrder(i)-1) + 1,6*(BeamSeg(segi).NodeOrder(i)-1)+1));
    end
end
xm=xm'/sm;

% Offset xm
% xm=xm-[0;0.2;0];

% Modify RB modes to comform to axes
for i = 1:NumNode
    Mode(Nmodes-5).Phi0(1:6,i) = [1;0;0;0;0;0];
    Mode(Nmodes-4).Phi0(1:6,i) = [0;1;0;0;0;0];
    Mode(Nmodes-3).Phi0(1:6,i) = [0;0;1;0;0;0];    
    Mode(Nmodes-2).Phi0(1:6,i) = [0;0;0;1;0;0];
    Mode(Nmodes-1).Phi0(1:6,i) = [0;0;0;0;1;0];
    Mode(Nmodes).Phi0(1:6,i)   = [0;0;0;0;0;1];
end

for i = 1:Nseg
    for j = 1:BeamSeg(i).NodeNo
        % velocity vector = angular velocity vector x displacement
        % Rotation produces displacement
        % Rotation in x axes --> no displacement
        Mode(Nmodes-2).Phi0(1:3,BeamSeg(i).NodeOrder(j)) = cross([1;0;0],BeamSeg(i).NodeX(j,:)'-xm);
        % Rotation in y axis --> displacement in z direction
        Mode(Nmodes-1).Phi0(1:3,BeamSeg(i).NodeOrder(j)) = cross([0;1;0],BeamSeg(i).NodeX(j,:)'-xm);
        % Rotation in z axis --> displacement in y direction
        Mode(Nmodes).Phi0(1:3,BeamSeg(i).NodeOrder(j))   = cross([0;0;1],BeamSeg(i).NodeX(j,:)'-xm);
    end
end

% Rewrite information back
for i = Nmodes-5:Nmodes
    for j = 1:6
        Vreal(j:6:end,i) = Mode(i).Phi0(j,1:NumNode)';
    end
end
%=========================================================================%

for i=Nmodes-5:Nmodes
    % ith mode
    k = i;
    for j = 1:6
        Mode(i).Omega = 0;  % (in rad/s), rigid body mode 
        Mode(i).Phi0(j,1:NumNode) = Vreal(j:6:end,i)'; % reshape storing in Phi0
        Mode(i).Phi1 = Mode(i).Phi0; % Velocity mode = displacement mode*w,  but w=0. Thus these are mass normalised
    end
    
    % This defines global values. Split-nodal values to be processed later.
    % RB velocity modes have same scaling as displacement modes, thus
    % multiply Vreal directly
    % Refer to Eqn. 5.19
    Mode(i).MPhi1 = Mmat*Vreal(:,i);
    Mode(i).MPhi1 = reshape(Mode(i).MPhi1,6,NumNode);

    % There is no Phi2 for RB modes!
  
  
    % Process split-nodal Phi0 and Phi1
    % Split end nodes, put other nodes into their respective places
    % DO NOT DIVIDE BY NUMBER OF CONNECTIONS
    % Also redefine coordinate directions
    Mode(k).Phi0S(1:6,NumIndices) = 0; % initialise
    Mode(k).Phi1S(1:6,NumIndices) = 0;
    for segi = 1:Nseg % Cycle through segments
      for j = 1:BeamSeg(segi).NodeNo % Cycle through nodes
          % current node
          Mode(k).Phi0S(:,BeamSeg(segi).NodeStart-1+j) = [BeamSeg(segi).GlobalAxes',zeros(3,3);zeros(3,3),BeamSeg(segi).GlobalAxes']*Mode(k).Phi0(:,BeamSeg(segi).NodeOrder(j));
          Mode(k).Phi1S(:,BeamSeg(segi).NodeStart-1+j) = [BeamSeg(segi).GlobalAxes',zeros(3,3);zeros(3,3),BeamSeg(segi).GlobalAxes']*Mode(k).Phi1(:,BeamSeg(segi).NodeOrder(j));   
      end
    end
end



%=========================================================================%
%                     Normal Structure Mode Shapes                        %
%=========================================================================%
for i = 1:Nmodes-6
    % ith mode
    k = i;
    for j = 1:6
        Mode(i).Omega = sqrt(Dreal(i));  % (in rad/s), extract from diagonal eigenvalue matrix    
        Mode(i).Phi0(j,1:NumNode) = Vreal(j:6:end,i)'; % take 1 of the 6 dof from each node, for ith mode
        Mode(i).Phi1 = Mode(i).Omega * Mode(i).Phi0; % Velocity mode = displacement mode * w (Eqn. 5.16)
    end
    % This defines global values. Split-nodal values to be processed later.
  
  
    % This MPhi1 defined in GLOBAL axes
    % Assume delta function shape function
    % Remember VReal is Phi0, therefore Phi1 is omegax
    Mode(i).MPhi1 = Mode(i).Omega*Mmat*Vreal(:,i);
    Mode(i).MPhi1 = reshape(Mode(i).MPhi1,6,NumNode);

    % Nodal force in global frame (Refer to Eqn 5.9)
    NodalForce = -Mmat*Vreal(:,i)*Dreal(i); % K*Vi(displacement eigenvector) = Force


    Mode(k).Phi2(1:6,1:NumIndices) = 0; % initialise force modes
    for segi = 1:Nseg % Cycle through segments
        for j = 1:BeamSeg(segi).NodeNo-1 % Cycle through nodes
            % Normal stress


            % find all child segments first
            % find according to depth level
            fg = 1;
            curseglist = [segi];
            totalchildseglist = [];          
            while fg > 0
                nextseglist = [];
                for nss = 1:numel(curseglist)
                    totalchildseglist = [totalchildseglist;BeamSeg(curseglist(nss)).ConnectTo];
                    nextseglist = [nextseglist;BeamSeg(curseglist(nss)).ConnectTo];
                end
                curseglist = nextseglist;
                if numel(curseglist) == 0
                    fg = 0;
                end
            end
            % totalchildseglist now contains all child segment number

            %=========================== Force ===========================%
            % plus all subsequent nodes on current segment
            for jj = j+1:BeamSeg(segi).NodeNo
                Mode(k).Phi2(1:3,BeamSeg(segi).NodeStart-1+j) = Mode(k).Phi2(1:3,BeamSeg(segi).NodeStart-1+j)+BeamSeg(segi).GlobalAxes'*NodalForce((BeamSeg(segi).NodeOrder(jj)-1)*6+1:(BeamSeg(segi).NodeOrder(jj)-1)*6+3);
            end

            % plus all child segments
            for kk = 1:numel(totalchildseglist) % For all child segments
                for ll = 2:BeamSeg(totalchildseglist(kk)).NodeNo % exclude first node on tree because it will be shared with parent
                    Mode(k).Phi2(1:3,BeamSeg(segi).NodeStart-1+j) = Mode(k).Phi2(1:3,BeamSeg(segi).NodeStart-1+j)+BeamSeg(segi).GlobalAxes'*NodalForce((BeamSeg(totalchildseglist(kk)).NodeOrder(ll)-1)*6+1:(BeamSeg(totalchildseglist(kk)).NodeOrder(ll)-1)*6+3);
                end
            end
            % This results in 1 less than the number of nodes. Defined at the
            % midpoint of each element

            %=============================================================%

            %========================== Moments ==========================%

            % plus all subsequent nodes on current segment
            for jj = j+1:BeamSeg(segi).NodeNo
                Mode(k).Phi2(4:6,BeamSeg(segi).NodeStart-1+j) = Mode(k).Phi2(4:6,BeamSeg(segi).NodeStart-1+j)+BeamSeg(segi).GlobalAxes'*NodalForce((BeamSeg(segi).NodeOrder(jj)-1)*6+4:(BeamSeg(segi).NodeOrder(jj)-1)*6+6);
                % and moments caused by forces
                % relative position vector x force vector
                Mode(k).Phi2(4:6,BeamSeg(segi).NodeStart-1+j) = Mode(k).Phi2(4:6,BeamSeg(segi).NodeStart-1+j)+BeamSeg(segi).GlobalAxes'*cross((BeamSeg(segi).NodeX(jj,:)'-(BeamSeg(segi).NodeX(j,:)'+BeamSeg(segi).NodeX(j+1,:)')/2),NodalForce((BeamSeg(segi).NodeOrder(jj)-1)*6+1:(BeamSeg(segi).NodeOrder(jj)-1)*6+3));
            end

            % plus all child segments
            for kk = 1:numel(totalchildseglist) % For all child segments
                for ll = 2:BeamSeg(totalchildseglist(kk)).NodeNo % exclude first node on tree because it will be shared with parent
                    Mode(k).Phi2(4:6,BeamSeg(segi).NodeStart-1+j) = Mode(k).Phi2(4:6,BeamSeg(segi).NodeStart-1+j)+BeamSeg(segi).GlobalAxes'*NodalForce((BeamSeg(totalchildseglist(kk)).NodeOrder(ll)-1)*6+4:(BeamSeg(totalchildseglist(kk)).NodeOrder(ll)-1)*6+6);
                    % and moments caused by forces
                    Mode(k).Phi2(4:6,BeamSeg(segi).NodeStart-1+j) = Mode(k).Phi2(4:6,BeamSeg(segi).NodeStart-1+j)+BeamSeg(segi).GlobalAxes'*cross((BeamSeg(totalchildseglist(kk)).NodeX(ll,:)'-(BeamSeg(segi).NodeX(j,:)'+BeamSeg(segi).NodeX(j+1,:)')/2),NodalForce((BeamSeg(totalchildseglist(kk)).NodeOrder(ll)-1)*6+1:(BeamSeg(totalchildseglist(kk)).NodeOrder(ll)-1)*6+3));
                end
            end
            %=============================================================%
        
        end
        
    end
  
    % Process split-nodal Phi0 and Phi1
    % Split end nodes, put other nodes into their respective places
    % DO NOT DIVIDE BY NUMBER OF CONNECTIONS
    % Also redefine coordinate directions
    Mode(k).Phi0S(1:6,NumIndices)=0; % initialise
    Mode(k).Phi1S(1:6,NumIndices)=0;
    for segi = 1:Nseg % Cycle through segments
      for j = 1:BeamSeg(segi).NodeNo % Cycle through nodes
          % Since Phi0 and Phi1 are in global coordinates, In beam's local coordinates, i.e. multiply by AxisMatrix'
          % current node
          Mode(k).Phi0S(:,BeamSeg(segi).NodeStart-1+j) = [BeamSeg(segi).GlobalAxes',zeros(3,3);zeros(3,3),BeamSeg(segi).GlobalAxes']*Mode(k).Phi0(:,BeamSeg(segi).NodeOrder(j));
          Mode(k).Phi1S(:,BeamSeg(segi).NodeStart-1+j) = [BeamSeg(segi).GlobalAxes',zeros(3,3);zeros(3,3),BeamSeg(segi).GlobalAxes']*Mode(k).Phi1(:,BeamSeg(segi).NodeOrder(j));

      end
    end
end




%=========================================================================%
%                     Obtain Curvatures Directly.                         %
%=========================================================================%
EMAT=[[     0     0     0      0      0       0     ];...
      [     0     0     0      0      0       0     ];...
      [     0     0     0      0      0       0     ];...
      [     0     0     0      0      0       0     ];...
      [     0     0    -1      0      0       0     ];...
      [     0     1     0      0      0       0     ]];


% C*phi2, which is strains and curvatures
% RB modes does not have phi2

%========================= Normal Structure Mode Shape ===================%
for k = 1:Nmodes-6
    % cycle through segments
    for segi = 1:Nseg
        L = BeamSeg(segi).NodeL(end);
        Xcoord = BeamSeg(segi).NodeL;
        % dL=L/NumInteg;
        % Xinteg=dL/2:dL:L-dL/2;

        Phi0k = zeros(6,BeamSeg(segi).NodeNo);
        Phi1k = zeros(6,BeamSeg(segi).NodeNo);
        Phi2k = zeros(6,BeamSeg(segi).NodeNo);
        
        % Select the part of Phi0S corresponding to the segment
        for i = 1:6
            for j = 1:BeamSeg(segi).NodeNo            
                Phi0k(i,j) = Mode(k).Phi0S(i,BeamSeg(segi).NodeStart-1+j);
                Phi1k(i,j) = Mode(k).Phi1S(i,BeamSeg(segi).NodeStart-1+j);
                Phi2k(i,j) = Mode(k).Phi2(i,BeamSeg(segi).NodeStart-1+j);
            end
        end

        % Multiply by nodelength for collocated cphi2. Refer to Eqn 5.21
        for j = 1:BeamSeg(segi).NodeNo-1  
            Mode(k).Seg(segi).cphi2(1:6,j) = (-((Phi0k(:,j+1)-Phi0k(:,j))/BeamSeg(segi).NodeDL(j) ...
                           -EMAT'*(Phi0k(:,j+1)+Phi0k(:,j))/2))*BeamSeg(segi).NodeDL(j);


        end
        
        % Put phi2 in this format too
        for i = 1:6
            Mode(k).Seg(segi).phi2 = Phi2k;
        end
        
        % Put phi1 in this format too
        % centrepoint, averaged
        for j = 1:BeamSeg(segi).NodeNo-1  
            Mode(k).Seg(segi).phi1(:,1:BeamSeg(segi).NodeNo-1) = (Phi1k(:,1:BeamSeg(segi).NodeNo-1)+Phi1k(:,2:BeamSeg(segi).NodeNo))/2;
        end
        
        
        %***********************ATTENTION*********************************%
        % Now Phi0, Phi1 and MPhi1 defined in global
        % Phi0S, Phi1S, Phi2 and cphi2 defined in local
        %*****************************************************************%
        
   end
end
%=========================================================================%


%=========================== Rigid Body Mode Shape =======================%
for k = Nmodes-5:Nmodes
    % cycle through segments
    for segi = 1:Nseg
        L = BeamSeg(segi).NodeL(end);
        Xcoord = BeamSeg(segi).NodeL;
        % dL=L/NumInteg;
        % Xinteg=dL/2:dL:L-dL/2;
        
        Phi0k = zeros(6,BeamSeg(segi).NodeNo);
        Phi1k = zeros(6,BeamSeg(segi).NodeNo);
        Phi2k = zeros(6,BeamSeg(segi).NodeNo);
        
        % Select the part of Phi0S corresponding to the segment
        for i = 1:6
            for j = 1:BeamSeg(segi).NodeNo
                Phi1k(i,j) = Mode(k).Phi1S(i,BeamSeg(segi).NodeStart-1+j);
            end
        end
        
        % Put phi1 in this format too
        % centrepoint, averaged
        for j = 1:BeamSeg(segi).NodeNo-1
            Mode(k).Seg(segi).phi1(:,1:BeamSeg(segi).NodeNo-1) = (Phi1k(:,1:BeamSeg(segi).NodeNo-1)+Phi1k(:,2:BeamSeg(segi).NodeNo))/2;
        end
    end
end
%=========================================================================%
% Now phi1, phi2, MPhi1, cphi2 are obtained.

%=========================================================================%
%            Compute coefficients of the nonlinear equations              %
%=========================================================================%

%****************** Attention, how would the modes be ordered ************%
%          Nmodes now refer to number of structural modes (minus RB)      %
%*************************************************************************%
Nmodes = Nmodes-6;
alpha1 = zeros(Nmodes+6,Nmodes+6);
alpha2 = zeros(Nmodes,Nmodes);

for j = 1:Nmodes+6
    k = j;
    alpha1(j,k) = 0;
    if j <= Nmodes
        if k <= Nmodes
            alpha2(j,k) = 0;
        end
    end
    % x1Mx1
    % In discrete form
    for N = 1:NumNode
        alpha1(j,k) = alpha1(j,k)+Mode(j).Phi1(:,N)'*Mode(k).MPhi1(:,N);
    end
    for segi = 1:Nseg
        for N = 1:BeamSeg(segi).NodeNo-1
            if j <= Nmodes && k <= Nmodes
                % x2Cx2
                alpha2(j,k) = alpha2(j,k)+Mode(j).Seg(segi).phi2(:,N)'*Mode(k).Seg(segi).cphi2(:,N);
            end
        end
    end
    % [j k]
    % end
end

% Scale modal coeffs
for j = 1:Nmodes
    mscale = sqrt(alpha1(j,j));
    cscale = sqrt(alpha2(j,j));
    for segi = 1:Nseg
        Mode(j).Seg(segi).phi1  = Mode(j).Seg(segi).phi1/mscale;
        Mode(j).Seg(segi).phi2  = Mode(j).Seg(segi).phi2/cscale;
        Mode(j).Seg(segi).cphi2 = Mode(j).Seg(segi).cphi2/cscale;      
    end
    Mode(j).Phi0  = Mode(j).Phi0/mscale;
    Mode(j).Phi1  = Mode(j).Phi1/mscale;
    Mode(j).MPhi1 = Mode(j).MPhi1/mscale;    
    Mode(j).Phi0S = Mode(j).Phi0S/mscale;
    Mode(j).Phi1S = Mode(j).Phi1S/mscale; 
    Mode(j).Phi2  = Mode(j).Phi2/cscale;      
end

for j = Nmodes+1:Nmodes+6
    mscale = sqrt(alpha1(j,j));
    for segi = 1:Nseg
        Mode(j).Seg(segi).phi1 = Mode(j).Seg(segi).phi1/mscale;      
    end
    Mode(j).Phi0  = Mode(j).Phi0/mscale;
    Mode(j).Phi1  = Mode(j).Phi1/mscale;
    Mode(j).MPhi1 = Mode(j).MPhi1/mscale;      
    Mode(j).Phi0S = Mode(j).Phi0S/mscale;
    Mode(j).Phi1S = Mode(j).Phi1S/mscale;       
end


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


% end

beta1  = zeros(Nmodes+6,Nmodes);
beta2  = zeros(Nmodes,Nmodes+6);
gamma1 = zeros(Nmodes+6,Nmodes+6,Nmodes+6);
gamma2 = zeros(Nmodes+6,Nmodes,Nmodes);

% Pre-process L1
for N = 1:NumNode
    for k = 1:Nmodes+6
        Mode(k).L1P1(:,:,N) = L1fun(Mode(k).Phi1(:,N)');
    end
end

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

% Gamma2 needed for deltacoeff
for segi = 1:Nseg
    for N = 1:BeamSeg(segi).NodeNo-1  
        for k = 1:Nmodes
            Mode(k).L2P2(:,:,segi,N) = L2fun(Mode(k).Seg(segi).phi2(:,N));
        end
    end
end

for segi = 1:Nseg
    for N = 1:BeamSeg(segi).NodeNo-1
        gm2 = zeros(Nmodes+6,Nmodes,Nmodes);
        
        for k = 1:Nmodes
            for l = 1:Nmodes
                L2kx2l = Mode(k).L2P2(:,:,segi,N) * Mode(l).Seg(segi).cphi2(:,N);
                for j = 1:Nmodes+6
                    % x1x2'+Ex2
                    % By definition, x2'+Ex2 = +-Mx1*w
                    % We know AinvB so do NOT compute beta
                    % Beta1 ant Beta2 are going to be defined via
                    % eigenvalues
                    
                    % x1L2(x2)Cx2
                    gm2(j,k,l) = gm2(j,k,l) + Mode(j).Seg(segi).phi1(:,N)'*L2kx2l;
                    
                end
            end
        end
        gamma2 = gamma2+gm2;
    end
end

% % Simple definition of alpha
% alpha1 = eye(Nmodes+6);
% alpha2 = eye(Nmodes);
% NO! Alpha is not diagonal any more
alpha1 = alpha1new;
alpha2 = alpha2new;

% beta1+
% beta2-
% q1. = wq2
% q2. = -wq1
EValue = zeros(Nmodes,1);

for i = 1:Nmodes
    EValue(i) = Mode(i).Omega;
end

beta1(1:Nmodes,:) = diag(EValue);
beta2(:,1:Nmodes) = -diag(EValue);

AMat = [alpha1,zeros(Nmodes+6,Nmodes);zeros(Nmodes,Nmodes+6),alpha2];

% BMat is W matrix in Eqn. 3.29
% Remember Eqn. 3.32
% alpha1_inverse * lambda1 = W_D
% alpha2_inverse * lambda2 = -W_D
BMat = [zeros(Nmodes+6,Nmodes+6),beta1;beta2,zeros(Nmodes,Nmodes)];


%=========================================================================%
%                Commence computation of AE coefficients                  %
%=========================================================================%
% Define global flight direction(in this case is the positive y axis)
VFGlobal = [0;1;0];

% Define Avector and bvector. Refer to eqn 2.122
% unsteady aerodynamic model, don't touch
nae = 2;
% Coefficients of a rational function approximation to Wagner's function
aeAvector = [0.165;0.335];
aebvector = [0.0455;0.3];
aeAb = aeAvector.*aebvector;

% define ae omega
aeomega = zeros(nae,1);
for i = 1:nae
    % aeomega(i)=aebvector(i)*Vi/aebvalue;
    % !!!A Vi is missing in here!!!
    aeomega(i)=aebvector(i)/aebvalue;
end

%Vector K_AE (Eqn 2.122b)
Vkt=[0,0,-1/aebvalue,(1-aeavalue),0,0];

% ==========This part need modification too!===========
% Compute transformation matrix from rb to velocity
% Must compute at c.g.!
btoV = zeros(6,6);
for i=Nmodes+1:Nmodes+6
    btoV(1:3,i-Nmodes) = Mode(i).Phi1(1:3,1);
    btoV(4:6,i-Nmodes) = Mode(i).Phi1(4:6,1);
end
btoVi = inv(btoV);
btoWonly=btoV(4:6,:);

% Compute rb flight speed
VComp=(btoV'*[VFGlobal;zeros(3,1)])';

% Matrix related to flap
muone   = zeros(Nmodes+6,Nmodes+6,Nmodes+6);
mutwo   = zeros(Nmodes+6,Nmodes+6,Nmodes+6);
muthree = zeros(Nmodes+6,Nmodes+6,Nmodes+6,Nflaps);

% Compute local aerodynamic orientation matrix
localCpM = zeros(6,6,Nseg);
for segi = 1:Nseg
    % Front and back
    % Define, y pointing along LE, z being normal and x
    % being the other direction
    % Transform into local axes before computing
    % coefficients (which is nondirectional)

    % Define local normal(z) direction. Assuming aerofoil
    % completely parallel to flow (!!!!!)
    % the local z direction is defined as the cross product
    % between flow direction and beam direction
    % The local x is then the cross product between y and z
    localY = VFGlobal;
    localX = BeamSeg(segi).GlobalAxes(:,1);
    localZ = cross(localX,localY);
    localZ = localZ/norm(localZ);
    localX = cross(localY,localZ);
    localX = localX/norm(localX);
    localCpM(:,:,segi) = [[localX,localY,localZ]',zeros(3,3);zeros(3,3),[localX,localY,localZ]'];
end


% Takes the local velocity at each node
A21 = A2fun_Arri(1, aeavalue, aebvalue);
for segi = 1:Nseg
    curseg = segi;
    % semi-DL
    NodeDLh = ([BeamSeg(segi).NodeDL/2;0]+[0;BeamSeg(segi).NodeDL/2]);
    % ******ADHOC REMOVAL OF CM
    % %     if segi<5
    % %         aecd0=0.01;
    % %         aecm0=0.025;
    % %     else
    % %         aecd0=0.02*2;
    % %         aecm0=0;
    % %     end
    localCp = localCpM(:,:,segi);
    
    CpPhi1   = zeros(6,Nmodes+6,BeamSeg(segi).NodeNo);
    A1CpPhi1 = zeros(6,6,Nmodes+6,BeamSeg(segi).NodeNo);
    A3CpPhi1 = zeros(6,6,Nmodes+6,BeamSeg(segi).NodeNo,Nflaps);
    
    if BeamSeg(segi).IsWing > 0
        % Take local aerodynamic values
        CL2 = BeamSeg(segi).CLA/2;
        CD0 = BeamSeg(segi).CD0;
        CM0 = BeamSeg(segi).CM0;
        % ===== use values provided here
        %         CLD=BeamSeg(segi).CLD;
        %         CMD=BeamSeg(segi).CMD;
        % pre-compute transformed velocities and matrices
        for k = 1:Nmodes+6
            for b = 1:BeamSeg(segi).NodeNo
                CpP = localCp*Mode(k).Phi1(:,BeamSeg(segi).NodeOrder(b));
                CpPhi1(:,k,b) = CpP;
                A1CpPhi1(:,:,k,b) = A1funCL(CpP,CL2,CD0,CM0,aeavalue,aebvalue);
                for nf = 1:Nflaps
                    A3CpPhi1(:,:,k,b,nf) = A3funD(CpP,CLD(nf),CMD(nf),aeavalue,aebvalue);
                end
            end
        end
        
        for i = 1:Nmodes+6
            for j = 1:Nmodes+6
                mmm = zeros(6,BeamSeg(segi).NodeNo);
                mcc = zeros(6,BeamSeg(segi).NodeNo);
                mdd = zeros(6,BeamSeg(segi).NodeNo,Nflaps);
                % pre-compute second-order couplings
                for b = 1:BeamSeg(segi).NodeNo
                    ndlhb = NodeDLh(b);
                    mmm(:,b) = (A1CpPhi1(:,:,i,b)*(CpPhi1(:,j,b)))*ndlhb;
                    % multiple flaps
                    for nf = 1:Nflaps
                        mdd(:,b,nf) = (A3CpPhi1(:,:,i,b,nf)*(CpPhi1(:,j,b)))*ndlhb;
                    end
                    % mcc(:,:,b)=(A21*CpPhi1(:,i,b)*(Vkt*CpPhi1(:,j,b)))*ndlhb;
                    % multiply by CL/2 (pi) here
                    mcc(:,b) = (A21*CpPhi1(:,i,b)*(Vkt*CpPhi1(:,j,b)))*ndlhb*CL2;
                end
                
                for k = 1:Nmodes+6
                    magm = 0;
                    magc = 0;
                    magd = zeros(Nflaps,1);
                    
                    
                    % Here, for now it is assumed that all lifting surfaces
                    % have the same design and size (which must be wrong)
                    
                    % Here, magc (mutwo) varies with airspeed!!! Different
                    % from the original paper where airspeed is included,
                    % thus here magc is a third order coupling term (Vi!)
                    
                    for b = 1:BeamSeg(segi).NodeNo
                        cpp1p = CpPhi1(:,k,b)';
                        magm = magm+cpp1p*mmm(:,b);
                        magc = magc+cpp1p*mcc(:,b);
                        for nf = 1:Nflaps
                            % <-----------------Flap Location Included in
                            % Here
                            magd(nf) = magd(nf)+cpp1p*mdd(:,b,nf)*flapdistr(nf,b);
                        end
                    end
                    muone(k,i,j) = muone(k,i,j)+magm;
                    mutwo(k,i,j) = mutwo(k,i,j)+magc;
                    
                    for nf = 1:Nflaps
                        muthree(k,i,j,nf) = muthree(k,i,j,nf)+magd(nf);
                    end
                    
                    % [segi,i,j,k]
                end
            end
        end
    end
end

% Scale by Rhoinf
% muone and mutwo's pi (CLA/2) multiplication is inside A1fun 
muone = muone*rhoinf*aebvalue;
mutwo = mutwo*2*rhoinf*aebvalue;            
muthree = muthree*rhoinf*aebvalue; 

%% Re-ordering nonlinear matrix coefficient

dcoeff = zeros(Nmodes,Nmodes+6,Nmodes);
for l = 1:Nmodes
    for j = 1:Nmodes
        for k = 1:Nmodes+6
            dcoeff(j,k,l) = gamma2(k,j,l);
        end
    end
end

goneord = zeros(Nmodes+6,Nmodes+6,Nmodes+6);
gtwoord = zeros(Nmodes,Nmodes,Nmodes+6);
gthreeord = zeros(Nmodes+6,Nmodes,Nmodes); %Delta 2 matrix coefficient(reshape of gtwoord)

moneord = zeros(Nmodes+6,Nmodes+6,Nmodes+6);
mtwoord = zeros(Nmodes+6,Nmodes+6,Nmodes+6);
mthreeord = zeros(Nmodes+6,Nmodes+6,Nmodes+6,Nflaps);

% reorder triples for optimisation
for k = 1:Nmodes+6
    for i = 1:Nmodes+6
        for j = 1:Nmodes+6
            goneord(i,j,k) = gamma1(k,i,j);
            moneord(i,j,k) = muone(k,i,j);
            mtwoord(i,j,k) = mutwo(k,i,j);
            for nf = 1:Nflaps
                mthreeord(i,j,k,nf) = muthree(k,i,j,nf);
            end
        end
    end
    for i = 1:Nmodes
        for j = 1:Nmodes
            gtwoord(i,j,k) = gamma2(k,i,j);
        end
    end
end

for k = 1:Nmodes
    for i = 1:Nmodes+6
        for j = 1:Nmodes
            gthreeord(i,j,k) = dcoeff(k,i,j);
        end
    end
end

StructMat.VrealNew = Vreal;

MC.alpha1 = alpha1;
MC.alpha2 = alpha2;
MC.beta1 = beta1;
MC.beta2 = beta2;
MC.goneord = goneord;
MC.gtwoord = gtwoord;
MC.gthreeord = gthreeord;
MC.moneord = moneord;
MC.mtwoord = mtwoord;
MC.mthreeord = mthreeord;
MC.AMat = AMat;
MC.BMat = BMat;

AeDef.aeomega = aeomega;
AeDef.aeAb = aeAb;
AeDef.nae = nae;
