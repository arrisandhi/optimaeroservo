function td=tilde(t)
% Cross product operator definition. Eqn 2.12
td = [0 -t(3) t(2);t(3) 0 -t(1);-t(2) t(1) 0]; 
end
% eof              