function [history,searchdir] = runfmincon_MF02_TM_MMax(maxMomOri,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,ShapeMass,BaseRA)
%% THIS FUNCTION PERFORM OPTIMISATION PROCESS 

fprintf(1, 'Run optimisation process.\n');
 
% Set up shared variables with OUTFUN
history.x    = [];
history.fval = [];
searchdir    = [];

% Fmincon configuration
X0 = ShapeMass;             % Initial value
LB = [-95/30.48;5];         % Lower bound of design vector
UB = [95/30.48;100];        % Upper bound of design vector

% Optimisation cost
Fun = @(X)main_optimCost_MF02_TM_MMax(X,BeamDef,BaseRA);

% Nonlinear constraint 
NLConstraintFun = @(X)confun_MF02_TM_MMax(X,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,maxMomOri);

% Set up of optimisation options
options = optimoptions(@fmincon,'OutputFcn',@outfun,'UseParallel',true,'Display','iter','Algorithm','active-set');

% Linear constraint AX = b
Atemp = [BeamDef.NodeL(1:BeamDef.Ndivs), ones(BeamDef.Ndivs,1)];
A = [-Atemp; Atemp];                                        
b = [-5*ones(BeamDef.Ndivs,1);100*ones(BeamDef.Ndivs,1)];

X = fmincon(Fun,X0,A,b,[],[],LB,UB,NLConstraintFun,options);

 function stop = outfun(x,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval]
           history.x = [history.x, x]
           fprintf('%d\n',history.fval);
           fprintf('%d\n',history.x);
         % Concatenate current search direction with 
         % searchdir.
         case 'done'
         otherwise
     end
 end
end