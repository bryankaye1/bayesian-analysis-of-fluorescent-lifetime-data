function [ebinary,estring] = read_matout_errorcheck(error,errparam)
ebinary = 0;
estring = '';
if sum(error)>0
    if error(1) == 1
        estring = sprintf('%s \nMarg error in PR vector loop 1\n',estring); %Continue this for other spots
    end
    if error(2) == 1
        estring = sprintf('%s \nMarg error in w02 vector loop 1\n',estring); %Continue this for other spots
    end
    if error(3) == 1
        estring = sprintf('%s \nMarg error in w2 vector loop 1\n',estring); %Continue this for other spots
    end
    if error(4) == 1
        estring = sprintf('%s \nMarg error in w1 vector loop 1\n',estring); %Continue this for other spots
    end
    if error(5) == 1
        estring = sprintf('%s \nMarg error in PR vector loop 2\n',estring); %Continue this for other spots
    end
    if error(6) == 1
        estring = sprintf('%s \nMarg error in w02 vector loop 2\n',estring); %Continue this for other spots
    end
    if error(7) == 1
        estring = sprintf('%s \nMarg error in w2 vector loop 2\n',estring); %Continue this for other spots
    end
    if error(8) == 1
        estring = sprintf('%s \nMarg error in w1 vector loop 2\n',estring); %Continue this for other spots
    end
    ebinary = 1;
end

if sum(errparam)>0
    if sum(errparam(1:5))>0
        estring = sprintf('%s \nThere is an error in errparam loop 1\n',estring); %Continue this for other spots
    end
    if sum(errparam(6:10))>0
        estring = sprintf('%s \nThere is an error in errparam loop 2\n',estring); %Continue this for other spots
    end
    if errparam(11)==111
        estring = sprintf('%s \nNaNs in loglike loop 1\n',estring);
    end
    if errparam(11)==123
        estring = sprintf('%s \nDid not search whole likelihood space loop 1\n',estring);
    end
    if errparam(12)==1
        estring = sprintf('%s \nAll prob in like are 0, loop 1\n',estring);
    end    
    if errparam(13)==111
        estring = sprintf('%s \nNaNs in loglike loop 2\n',estring);
    end
    if errparam(13)==123
        estring = sprintf('%s \nDid not search whole likelihood space loop 2\n',estring);
    end
    if errparam(14)==1
        estring = sprintf('%s \nAll prob in like are 0, loop 2\n',estring);
    end     
    ebinary =1;  
end
if ebinary==1
    fprintf('\n%s\n',estring);
end

end