function A3=A3funD(x1, CLD, CMD, aeav, aebv)

A3=[[0,0,0,0,0,0]; ...
    [0,-CLD*x1(3),0,0,0,0]; ...
    [0,CLD*x1(2),0,0,0,0]; ...
    [0,CLD*x1(2)*aeav*aebv+CMD*x1(2)*aebv*2,0,0,0,0]; ...
    [0,0,0,0,0,0]; ...
    [0,0,0,0,0,0]];

                    