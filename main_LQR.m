function [u,x,y,Klqr] = main_LQR(ssLinRed,StateWORBRed,SimDef,GustDef)
%% PLANT MODEL
% Set discrete time MPC controller dynamics
SrA = ssLinRed.a;
SrB = ssLinRed.b;
SrC = ssLinRed.c;
SrD = ssLinRed.d;
dt = SimDef.dt;             % Time step
tSim = SimDef.tSim;         % Simulation time
[Ad,Bd,Cd,~]=c2dm(SrA,SrB,SrC,SrD,dt); % Continous to discrete system

BdG = Bd(:,end);        % Separate gust input from input matrix
Bd = Bd(:,1:end-1);

n = size(Ad,1);     % Ad row size 
m = size(Bd,2);     % Bd column size
p = size(Cd,1);     % Cd row size

%% SET UP CONTROL WEIGHT MATRIX
Q = 20*diag(ones(1,n));      % State weight matrix
R = 50*diag(ones(1,m));      % Input weight matrix  

%% SET UP DLQR
[Klqr,S,cl_poles]=dlqr(Ad,Bd,Q,R);

ul = -0.5*ones(m,1);        % Min deflection constraint
uh = 0.5*ones(m,1);         % Max deflection constraint

%% SIMULATION LOOP
x(:,1) = StateWORBRed;      
u = zeros(m,length(tSim));
y = zeros(p,length(tSim));

Vgust_alltime = GustDef.Vgust_alltime; 

for it = 1:numel(tSim); % Simulation
    %% DLQR
    u(:,it) = -Klqr*x(:,it);
    
    % Introduce control saturation
    if u(1,it)> uh(1)
       u(1,it) = uh(1);
    elseif u(1,it) < ul(1)
       u(1,it) = ul(1);
    end
 
    x(:,it+1) = Ad*x(:,it)+Bd*u(:,it)+BdG*Vgust_alltime(3,it);
    y(:,it) = Cd*x(:,it);  
end

    
    
 