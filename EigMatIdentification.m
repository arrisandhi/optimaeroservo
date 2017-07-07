% This code plots poles of EigMat 

% Initialise CurrentState
CurrentState = trimstate;

% Add eigenvector as state to corresponding CurrentState
for iEval = 1
    CurrentState(1:RVnum)                   = EigMatEvec(1:RVnum,iEval);
    CurrentState([1:RVnum]+neig+6)          = EigMatEvec(RVnum+[1:RVnum],iEval);
    CurrentState([1:RVnum]+2*neig+6)        = EigMatEvec(2*RVnum+[1:RVnum],iEval);
    CurrentState([1:RVnum]+2*neig+6+neig+6) = EigMatEvec(3*RVnum+[1:RVnum],iEval);
    CurrentState = real(CurrentState);
    [ModalRE1, ModalCE1, dRH1, dCH1, dR1, dC1] = dispintegrationvar_Arri(CurrentState, nInterest, ttime, NumNode, neig, NumIndices, BeamSeg, Nseg, Mode);
    h = figure('Visible','off');
    plot(ModalRE1(1,:),ModalRE1(3,:))
    saveas(h,sprintf('%d.jpg',iEval));
end


