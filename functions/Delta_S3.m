function D_S = Delta_S3(r,q,aff2)
    % q: (n,1) or (1,n) vector
    % r: (6,1) or (1,6) vector   [r11 r21 theta1 r12 r22 theta2] 
    % aff2: (9,1) vector, atomic form factor .^2
    % r(3) = sqrt(r(1)^2 + r(2)^2 -2*r(1)*r(2)*cos(r(3)));
    % r(6) = sqrt(r(4)^2 + r(5)^2 -2*r(4)*r(5)*cos(r(6)));

    D_S = r(7) * 2 * aff2 .*(math_sinc(q*r(4)) + math_sinc(q*r(5)) + math_sinc(q*r(6)) - math_sinc(q*r(1)) - math_sinc(q*r(2)) - math_sinc(q*r(3))); 


end

