function [ssLinRed,StateWORBRed,Ti,hsv] = MOR(wingSSG,StateWORB)

[ssLinBal,hsv,Tr,Ti] = balreal(wingSSG);  % Balanced realisation

% Calculate importance level of states via grammian vector evaluation
ImpLvl = hsv/sum(hsv);

iRed = find(ImpLvl < 0.001);

% Index of reduced states
ssLinRed = modred(ssLinBal,iRed,'Truncate');        % Model Reduction
StateWORBRed = Tr*StateWORB;                        % Transforms the states into reduced states
StateWORBRed = StateWORBRed(1:length(ssLinRed.a));

end