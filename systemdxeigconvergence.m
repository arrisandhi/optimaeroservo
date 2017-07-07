function EigMatFull = systemdxeigconvergence(MC,AeDef,trimstate,btoV,Nmodes)
% This function returns the temporary dynamic matrix (EigMat). This EigMat still needs to be
% multiplied by Ainverse from guyancomp function
%
% INPUT:
    % trimstate & 2*Nmodes+6+nae*(Nmodes+6) x 1 		& Trim state    \\ \hline
    % neig		& 1x1 	 	& Number of modes                           \\ \hline
    % nae		& 1x1      	& Number of aerodynamic state               \\ \hline
    % aeomega	& nae x 1   & aebvector divided by aebvalue             \\ \hline
    % aeAb	 	& nae x 1       & A multiply by b (Theodorsen solution) \\ \hline
    % goneord	& Nmodes+6 x Nmodes+6 x Nmodes+6 & Matrix gamma1        \\ \hline
    % gtwoord	& Nmodes x Nmodes x Nmodes+6     & Matrix gamma2 matrix \\ \hline
    % moneord	& Nmodes+6 x Nmodes+6 x Nmodes+6 & Matrix H1-flap       \\ \hline
    % mtwoord 	& Nmodes+6 x Nmodes+6 x Nmodes+6 & Matrix H2-flap       \\ \hline	
    % btoV		& 6x6                 & Transformation matrix from rb to velocity \\ \hline
    % Bmat      & Nmodes+6 x Nmodes+6 & Full Matrix WD \\ \hline

% OUTPUT:
    % EigMatFull& 1x1 & Dynamic matrix \\ \hline

nae = AeDef.nae;
aeomega = AeDef.aeomega;
aeAb = AeDef.aeAb;
goneord = MC.goneord;
gthreeord = MC.gthreeord;
moneord = MC.moneord;
mtwoord = MC.mtwoord;
BMat = MC.BMat;

% Velocity state and Force State
x = trimstate(1:2*Nmodes+6);

% Aerodynamic state
xae = trimstate(2*Nmodes+6+1:2*Nmodes+6+nae*(Nmodes+6));

% Convert ae states to matrix form for better manipulation
xam = zeros(Nmodes+6,nae*(Nmodes+6));
for i = 1:nae
    xam(:,(Nmodes+6)*(i-1)+1:(Nmodes+6)*i) = eye(Nmodes+6)*aeAb(i);
end
xaa = xam*xae;

% Velocity state
x1 = x(1:Nmodes+6);

% Initialisation
EigMatLin = zeros(2*Nmodes+6+nae*(Nmodes+6),2*Nmodes+6+nae*(Nmodes+6));
EigMatLin(1:2*Nmodes+6,1:2*Nmodes+6) = BMat; % Obtain Lambda matrix

EigMat = zeros(2*Nmodes+6+nae*(Nmodes+6),2*Nmodes+6+nae*(Nmodes+6));


% Structural nonlinear coeff, does not contribute as effectively x1=0
% populate matrix
for k = 1:Nmodes+6
    EigMat(k,1:Nmodes+6) = EigMat(k,1:Nmodes+6)-x1'*goneord(:,:,k)-x1'*goneord(:,:,k)';
end

for k = 1:Nmodes
    EigMat(Nmodes+6+k,Nmodes+7:2*Nmodes+6) = EigMat(Nmodes+6+k,Nmodes+7:2*Nmodes+6)+x1'*gthreeord(:,:,k);
end

% Rather, Compute current flight speed Vi (all dirs flight)
% Compute aerodynamic state influences
btv = btoV(1:3,:);          % Transformation matrix from rb to velocity
bsq = btv'*btv;            
Vi = sqrt(x1(Nmodes+1:end)'*bsq*x1(Nmodes+1:end));

% Refer to Eqn. 3.42
for k = 1:Nmodes+6
    % H1 on x1
    EigMat(k,1:Nmodes+6) = EigMat(k,1:Nmodes+6)+x1'*(moneord(:,:,k))+x1'*(moneord(:,:,k))';
    % H2 on x1xz
    EigMat(k,1:Nmodes+6) = EigMat(k,1:Nmodes+6)+(mtwoord(:,:,k)*(xaa*Vi))';
    % H2 on Vi
    EigMat(k,Nmodes+1:Nmodes+6) = EigMat(k,Nmodes+1:Nmodes+6)+(x1)'*(mtwoord(:,:,k)*xaa*x1(Nmodes+1:end)'*bsq/Vi);
    % H2 on aero
    EigMat(k,2*Nmodes+6+1:2*Nmodes+6+nae*(Nmodes+6)) = EigMat(k,2*Nmodes+6+1:2*Nmodes+6+nae*(Nmodes+6))+x1'*mtwoord(:,:,k)*Vi*xam;
end

% Finally, influence of Vinf on aero states
for k = 1:Nmodes+6
    for i = 1:nae
        % x1 component
        EigMatLin(2*Nmodes+6+(i-1)*(Nmodes+6)+k,k) = EigMatLin(2*Nmodes+6+(i-1)*(Nmodes+6)+k,k)+1;   
        % ae component
        EigMatLin(2*Nmodes+6+(i-1)*(Nmodes+6)+k,2*Nmodes+6+(i-1)*(Nmodes+6)+k) = EigMatLin(2*Nmodes+6+(i-1)*(Nmodes+6)+k,2*Nmodes+6+(i-1)*(Nmodes+6)+k)-aeomega(i)*Vi;
        % Vi component
        if k > Nmodes
            EigMatLin(2*Nmodes+6+(i-1)*(Nmodes+6)+k,Nmodes+1:Nmodes+6) = EigMatLin(2*Nmodes+6+(i-1)*(Nmodes+6)+k,Nmodes+1:Nmodes+6)-xae((i-1)*(Nmodes+6)+k)*aeomega(i)*x1(Nmodes+1:end)'*bsq/Vi;
        end
    end
end

EigMatFull = EigMatLin + EigMat;