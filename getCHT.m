function CHT = getCHT(Cht, epsilon)
    if strcmp(Cht, 'EC')
        CHT.str = Cht;
        CHT.epsilon = epsilon.epsilon0;
    else
        CHT.str = Cht;
        CHT.epsilon = 0;
    end
end