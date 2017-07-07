function [Cmp,Rvp] = RotDisp(tSpan,ModalCE,ModalRE,CurrentState)

Vstate = CurrentState(1:neig+6);

Cmp = zeros(3,3,NumNode,numel(tSpan));
Rvp = zeros(3,NumNode,numel(tSpan));
Cmp(:,:,:,1) = ModalCE;
Rvp(:,:,1) = ModalRE;

for t = 2:numel(tSpan)
    Vt = zeros(6,NumNode);
    for m = 1:neig+6
        Vt = Vt + Mode(m).Phi1*Vstate(m);
    end
    
    for i = 2:NumNode
        Cm = ModalCE(:,:,i);
        Cmp(:,:,i,t) = Cm*tilde(Vt(4:6,i));
        Rvp(:,i,t) = Cm*Vt(1:3,i);
    end
end
