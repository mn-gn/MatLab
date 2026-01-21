function S = Debye3(q,r,coeff) 
    % q: (n,1) or (1,n) vector
    % r: (3,1) or (1,3) vector   [r1 r2 theta] 
    % atom: string ex('H1-')


    f = aff(coeff,q);
    f2 = f.^2;

    R = zeros(3);
    R(find(triu(ones(3), 1))) = r;
    R = R + R';

    S = 0;
    for i = 1 : 3
      for j = 1 : 3
            S = S + f2 .* math_sinc(q*R(i,j));

       end
    end

end