function L2=L2fun(x2)

L2=zeros(6,6); 
F=zeros(3,1);
M=zeros(3,1);

F(1:3)=x2(1:3);
M(1:3)=x2(4:6);

L2(4:6,4:6)= [[  0   -M(3)  M(2)]; ...
              [ M(3)  0    -M(1)]; ...
              [-M(2)  M(1)  0   ]];

          
L2(1:3,4:6)= [[  0   -F(3)  F(2)]; ...
              [ F(3)  0    -F(1)]; ...
              [-F(2)  F(1)  0   ]];

             
L2(4:6,1:3)= L2(1:3,4:6);                          
% eof              