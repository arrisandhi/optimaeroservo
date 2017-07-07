% figure;
length = 0:BeamDef.Length/100:BeamDef.Length;
tSpan = 0:10/100:10;
tgi = 1;
tge = 1+2;
uref = 3;
Ly = 1;

for s = 1:numel(length)
    for it = 1:numel(tSpan)
        if tSpan(it) <= tge && tSpan(it) >= tgi
            Vgust(s,it) = uref/2*(1-cos(2*pi*(tSpan(it)-tgi)/(tge-tgi)))*1/2*cos(pi*(length(s))/(10)); % Gust function
%            Vgust(s,it) = -1/2*cos(pi*(length(s))/(10)); % Gust function
        else
            Vgust(s,it) = 0;
        end
    end
end

surf(length(1:end),tSpan(11:end),Vgust(1:end,11:end)')
axis equal
% xlabel('tSpan');ylabel('length');zlabel('vg')