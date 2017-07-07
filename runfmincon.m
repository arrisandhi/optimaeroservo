function [history,searchdir] = runfmincon(StrainEnergyOri,maxZDisp,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,MC,Mode,Nmodes,btoV,RVnum)
 fprintf(1, 'Run optimisation process.\n');
 
% Set up shared variables with OUTFUN
history.x    = [];
history.fval = [];
searchdir    = [];

% Fmincon configuration
X0 = BeamProp.RA;           % Initial value
LB = 5*ones(20,1);          % Lower bound of design vector
UB = 100*ones(20,1);        % Upper bound of design vector

Fun = @(RA)main_optimCost(RA,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,Nmodes,RVnum,LB,UB,maxZDisp);
NLConstraintFun = @(X)confun(X,BeamProp,BeamDef,BeamSeg,AeDef,GustDef,SimDef,MC,Mode,Nmodes,btoV,RVnum,maxZDisp,StrainEnergyOri);
options = optimoptions(@fmincon,'OutputFcn',@outfun,'UseParallel',true,'Display','iter','Algorithm','active-set');

X = fmincon(Fun,X0,[],[],[],[],LB,UB,NLConstraintFun,options);

 function stop = outfun(x,optimValues,state)
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval]
           history.x = [history.x, x]
           fprintf('%d\n',history.fval);
           fprintf('%d\n',history.x);
         % Concatenate current search direction with 
         % searchdir.
           searchdir = [searchdir;... 
                        optimValues.searchdirection'];
           plot(x(1),x(2),'o');
         % Label points with iteration number and add title.
         % Add .15 to x(1) to separate label from plotted 'o'
           text(x(1)+.15,x(2),... 
                num2str(optimValues.iteration));
           title('Sequence of Points Computed by fmincon');
         case 'done'
             hold off
         otherwise
     end
 end
end