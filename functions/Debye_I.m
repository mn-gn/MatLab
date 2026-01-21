function S = Debye_I(q,r)
% r = [r1, r2, theta]
n = length(r);

aff_coeffs_I = [20.1472 4.347 18.9949 0.3814 7.5138 27.766 2.2735 66.8776 4.0712];
aff_I = aff(aff_coeffs_I,q);
f2I= aff_I.^2;

r(3) = sqrt(r(1)^2 + r(2)^2 - 2*r(1)*r(2)*cos(r(3)));

R = zeros(n);
R(find(triu(ones(n), 1))) = r;
R = R + R';

S = 0;
for i = 1 : n
    for j = 1 : n
        S = S + f2I .* math_sinc(q*R(i,j));

    end
end

end

