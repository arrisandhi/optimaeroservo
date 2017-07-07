for i = 55:0.1:60
    AeDef.Vi = i;
    main_pre_Goland;
    main_processing_Goland;
    
    if max(real(eig(wingSSG.a))) > 0
        break
    else
        plot(eig(wingSSG.a),'o')
        hold on;
    end
    i
end
