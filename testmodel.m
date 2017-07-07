main_pre_Goland_02;
main_processing_Goland;
[u_20_5,du_20_5,x_20_5,y_20_5] = main_MPCsptwir(ssLinRed,StateWORBRed,SimDef,GustDef);

zDispWORBRed_20_5 = calc_disprot(u_20_5,x_20_5,y_20_5,BeamDef,SimDef,trimstate,StateWORBRed,StateWORB,Nmodes,BeamSeg,Mode,Ti,RVnum);