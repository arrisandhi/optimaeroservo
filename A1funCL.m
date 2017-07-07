function A1=A1funCL(x1, CL2, CD0, CM0, aeav, aebv)

% Input CL/2, CD0, CM0, aeavalue and aebvalue
A1=[[0,0,0,0,0,0]; ...
    [0,-CD0*x1(2),x1(3)*CL2,-aebv*(1-aeav)*x1(3)*CL2,0,0]; ...
    [0,-x1(3)*CL2,0,aebv*(1-aeav)*x1(2)*CL2,0,0]; ...
    [0,-aeav*aebv*x1(3)*CL2+aebv*2*CM0*x1(2),0,(aeav-aeav*aeav-0.5)*aebv*aebv*x1(2)*CL2,0,0]; ...
    [0,0,0,0,0,0]; ...
    [0,0,0,0,0,0]];