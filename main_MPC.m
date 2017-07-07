function [u,du,x,y] = main_MPC(ssLinRed,StateWORBRed,SimDef,GustDef)
% MPC set point tracking including input rate cost. This code performs tracking of piecewise constant
% reference signal(uref, xref). 
% Cost function is a function of :
% 1. the difference between the state value and the state reference signal,
% 2. the difference between the current input value and the reference input.
% 3. the difference between the next input and current input (input rate)

% The arrangement of quadratic programming matrix is as follow:
% G = diag([RSQ RSQ RSQ ... RSP])
% other QP matrix follow the arrangement of G

%% PLANT MODEL
% Set discrete time MPC controller dynamics
SrA = ssLinRed.a;
SrB = ssLinRed.b;
SrC = ssLinRed.c;
SrD = ssLinRed.d;
dt = SimDef.dt;                         % Time step
tSim = SimDef.tSim;                     % Simulation time
[Ad,Bd,Cd,~]=c2dm(SrA,SrB,SrC,SrD,dt);  % Continous to discrete system

BdG = Bd(:,end);        % Separate gust input from input matrix
Bd = Bd(:,1:end-1);

n = size(Ad,1);     % Ad row size 
m = size(Bd,2);     % Bd column size
p = size(Cd,1);     % Cd row size

%% SET UP QP
P = 20*diag(ones(1,n));      % Terminal weight matrix
Q = 20*diag(ones(1,n));      % State weight matrix
R = 50*diag(ones(1,m));      % Input weight matrix 
S = 50*diag(ones(1,m));      % Input rate weight matrix

N = 50; % Prediction horizon

% Assembly matrix for Quadratic Programming
RSQ = blkdiag(R,S,Q);
RSP = blkdiag(R,S,P);
G = [kron(eye(N-1),RSQ) zeros(size(RSQ,1)*(N-1),size(RSP,2)); zeros(size(RSP,1),(N-1)*size(RSQ,2)) RSP];

% General expression of QP is 1/2*x'*H*x+F'x. In this case F matrix is zero
H = 2*G;

uRef = zeros(m,1);      % Input reference 
xRef = StateWORBRed;    % State reference

%% CONSTRAINTS FOR QP
% Equality constraints
DI = [-Bd zeros(n,m) eye(n);...
      -eye(m) eye(m) zeros(m,n)];
DIr = size(DI,1);
DIc = size(DI,2);
DDIFull = [kron(eye(N-1),DI) zeros((N-1)*DIr,DIc);...
       zeros(n,N*DIc)];

OI = [zeros(n,m) zeros(n,m) zeros(n,n);...
      eye(m) zeros(m,m) zeros(m,n)];
DOI = kron(eye(N-1),OI);
DOIFull = [zeros((N-1)*DIr,DIc) DOI ;...
           zeros(n,N*DIc)];

OA = [zeros(n,m) zeros(n,m) -Ad;...
      zeros(m,m) zeros(m,m) zeros(m,n)];
DOA = kron(eye(N-2),OA);
DOAFull = [zeros(DIr,N*DIc);...
          DOA zeros((N-2)*DIr,2*DIc);...
          zeros(n,(N-2)*DIc) zeros(n,m) zeros(n,m) -Ad -Bd zeros(n,m) eye(n)];

D = DDIFull + DOIFull + DOAFull;


% Inequality constraints
E1 = [eye(m) zeros(m,m) zeros(m,n);...
     -eye(m) zeros(m,m) zeros(m,n);...
     zeros(m,m) eye(m) zeros(m,n) ;...
     zeros(m,m) -eye(m) zeros(m,n) ;...
     zeros(p,m) zeros(p,m) Cd ;...
     zeros(p,m) zeros(p,m) -Cd];
E = kron(eye(N),E1);

ul = -0.348*ones(m,1);      % Input lowerbound
uh = 0.348*ones(m,1);       % Output upperbound
dul = -4.36e-3*ones(m,1);   % Input rate lowerbound
duh = 4.36e-3*ones(m,1);    % Input rate upperbound
yl = -1000000*ones(p,1);    % Output lowerbound
yh = 1000000*ones(p,1);     % Output upperbound

h1 = [uh+uRef; -ul+uRef; duh; -dul; yh-Cd*xRef; -yl+Cd*xRef];
h = repmat(h1,N,1);


%% SIMULATION LOOP
options = optimoptions('quadprog','Display','iter');
x(:,1) = StateWORBRed;
u = zeros(m,length(tSim));
y = zeros(p,length(tSim));

Vgust_alltime = GustDef.Vgust_alltime;

for it = 1:numel(tSim) - N;
    
    gLIDAR = Ad*x(:,it) - xRef + Bd*uRef + BdG*Vgust_alltime(3,it); 
    for p = 1:N-1
        gLIDAR = [gLIDAR;...
            zeros(m,1);...
            (Ad-eye(n))*xRef + Bd*uRef + BdG*Vgust_alltime(3,it+p)];
    end
    
    % Solve QP problem using quadprog
    X = quadprog(H,[],E,h,D,gLIDAR,[],[],[],options);
    u(:,it)  = X(1:m); + uRef;
    du(:,it) = X(m+1:m+m);
    x(:,it+1)= X(2*m+1:2*m+n) + xRef;
    y(:,it)  = Cd*x(:,it);
end

    
    
 