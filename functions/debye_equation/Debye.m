function S = Debye(q,r,aff_c)
% q = scttering vector
% r = [r1, r2, r3]

n = length(r);      % n = 1, 3, 6, 10, ...
f = aff(aff_c,q);

f2= f.^2;

R = zeros(n);
R(find(triu(ones(n), 1))) = r;

R = R + R';

S = 0;
for i = 1 : n
    for j = 1 : n
        S = S + f2 .* math_sinc(q*R(i,j));

    end
end

end

