% THIS CODE CONTAINS ALL INPUTS OF PROGRAM

%% BEAM DEFINITION
BeamDef.Nbeam   = 1;
BeamDef.Ndivs   = 20;                               % Number of segments between each mass node
BeamDef.NumNode = BeamDef.Ndivs+1;                  % Number of nodes
BeamDef.Nseg    = 1;                                % For wing with multiple segment
BeamDef.NumIndices = BeamDef.NumNode+BeamDef.Nseg-1;% Number of actual indices, including duplications at connections
BeamDef.Length  = 30.48;                               % Length of beam (semi-span)
BeamDef.Ndl     = BeamDef.Length/BeamDef.Ndivs;     % Node-to-node distance everywhere
BeamDef.NodeL   = (0:BeamDef.Ndl:BeamDef.Ndivs*BeamDef.Ndl)';

%% CONSTRUCT STRUCTURE DATA TYPE OF BEAM SEGMENT
BeamSeg(1).EndConn   = 1;                                   %????????
BeamSeg(1).ConnectTo = [];                                  %????????
BeamSeg(1).Parent    = 0;                                   %????????
BeamSeg(1).NodeNo    = BeamDef.Ndivs+1;
BeamSeg(1).NodeOrder = (1:BeamDef.Ndivs+1)';
BeamSeg(1).NodeStart = 1;
BeamSeg(1).NodeDL    = ones(BeamDef.Ndivs,1)*BeamDef.Ndl;                   % Length of each segment
BeamSeg(1).NodeL     = (0:BeamDef.Ndl:BeamDef.Ndivs*BeamDef.Ndl)';                  % Coordinate of mass node along wing
BeamSeg(1).NodeX     = [BeamSeg(1).NodeL,zeros(BeamDef.Ndivs+1,2)]; % Coordinate of mass node along wing in 3 axis(xyz frame)
BeamSeg(1).GlobalAxes = eye(3);


%% MASS DISTRIBUTION SHAPE FUNCTION
ShapeMass = X0(1:3);

%% CALCULATE BEAM PROPERTIES 
BeamProp = main_beamprop_MF03CW(BeamDef,ShapeMass);

%% AERODYNAMIC CHARACTERISTIC OF WING
BeamSeg(1).IsWing   = 1;            
BeamSeg(1).CLA      = 2*pi;
BeamSeg(1).CD0      = 0.0;
BeamSeg(1).CM0      = 0.0;
% Initialisation of flap aerodynamics coefficient along wing
BeamSeg(1).CLD = [];
BeamSeg(1).CMD = [];

% Ratio of distance of aerodynamic center(a.c) from elastic axis(e.a) to 
% half chord. If e.a. is backwards of a.c. then a is positive. The distance is a*b
AeDef.aeavalue = 0.16;                
AeDef.aebvalue = 0.9144;              % half-chord b

%% FLAP AERODYNAMICS CHARACTERISTICS AND POSITION
AeDef.Nflaps = 1;
AeDef.CLD    = [1,1];
AeDef.CMD    = [-0.25,-0.25];
AeDef.trimFlap = 0;    % Trim position of flap
AeDef.flapdistr = zeros(AeDef.Nflaps,BeamSeg(1).NodeNo);
AeDef.flapdistr(1,13:21) = 1;                                                     

%% FREESTREAM FLOW PROPERTIES
AeDef.Vi       = 28;               % Free stream velocity
AeDef.rhoinf   = 1.02;             % Air density

%% Simulation time 
SimDef.tStart = 0;      % Time Start
SimDef.tEnd   = 15;     % Time End
SimDef.dt     = 0.05;   % Timestep
SimDef.tSim   = SimDef.tStart:SimDef.dt:SimDef.tEnd;

%% Select number of considered modes
% Number of retained modes which are used for calculating nonlinear matrix
% and constructing dynamic matrix of aeroservoelastic system
Nmodes = 20;

%% Cantilever wing state space parameter
% For the case of cantilever wing, the states which are related to rigid
% body motion must be excluded when constructing the dynamic matrix. RVnum
% is the number of state exclude rigid body state.
RVnum = 20;

%% DEFINE GUST DISTRIBUTION PARAMETER
GustDef.Ly   = 1;         % Spatial distribution parameter
GustDef.uref = 17.07;     % Gust velocity reference(max value)
GustDef.tgi  = 3;         % Gust starting time
GustDef.tge  = 3+0.7143 ;     % Gust ending time

%% CALCULATE VGUST AND INPUT VECTOR ALONG SIMULATION TIME
GustDef.Vgust_alltime = zeros(3,numel(SimDef.tSim));
for it = 1:numel(SimDef.tSim)
    if SimDef.tSim(it) <= GustDef.tge && SimDef.tSim(it) >= GustDef.tgi
        GustDef.Vgust_alltime(3,it) = GustDef.uref/2*(1-cos(2*pi*(SimDef.tSim(it)-GustDef.tgi)/(GustDef.tge-GustDef.tgi))); % Gust function
    else
        GustDef.Vgust_alltime(3,it) = 0;
    end
    
end
