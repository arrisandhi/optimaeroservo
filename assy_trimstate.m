function trimstate = assy_trimstate(AeDef,Nmodes,btoV)
%% Assembly trimstate
% INPUTS:
    % nae     : Number of aerodynamic states
    % aeomega : aebvector/aebvalue
    % Vi      : Free stream velocity
  
% OUTPUTS. The state of system is assembled as follow:
    % trimstate          1 to Nmodes     : velocity(velocity and angular velocity) state structural
    % trimstate   Nmodes+1 to Nmodes+6   : velocity(velocity and angular velocity) state rigid body
    % trimstate   Nmodes+7 to 2*Nmodes+6 : force moment state
    % trimstate 2*Nmodes+7 to  end       : aerodynamics state

nae = AeDef.nae;
aeomega = AeDef.aeomega;
Vi = AeDef.Vi;

% associated with each velocity mode
trimstate = zeros(2*Nmodes+6+nae*(Nmodes+6),1);

% trimsate(Nmodes+2,1) : velocity mode rigid body in y direction 
% btoV(2,2) : velocity mode in y direction
trimstate(Nmodes+2,1) = Vi/btoV(2,2);

% Definition of velocity state
x1 = trimstate(1:Nmodes+6);

% Initialisation of aerodynamic state(written in nae column)
xaemat = zeros(Nmodes+6,nae);

% Write Back AE STATES
for k = 1:Nmodes+6
    for i = 1:nae
        xaemat(k,i) = x1(k)/Vi/aeomega(i,1);
    end
end

% Reshape aerodynamic state matrix and put them into trimstate vector
xae = reshape(xaemat,nae*(Nmodes+6),1);
trimstate(2*Nmodes+6+1:2*Nmodes+6+nae*(Nmodes+6)) = xae;

