%%%atomic form factor

function f = aff(pars,q)
    f = pars(9);

    for i = 1 : 4
        f = f + pars(2*i-1)*exp(-pars(2*i)*(q/(4*pi)).^2);

    end


end