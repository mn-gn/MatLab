S = zeros(size(S_data_all));

for i = 1 : 101
    r = [final_r123; final_locals(:,i); final_r7];
    S(:,i) = Delta_S3(r,q_data,f2);
    
end

asdf = final_locals;
for i = 1: 101

asdf(3,i) = sqrt(final_locals(1,i)^2 + final_locals(2,i)^2 - 2*final_locals(1,i)*final_locals(2,i)*cos(final_locals(3,i)));
asdf(:,i) = sort(asdf(:,i));

end