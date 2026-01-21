clear; close all;

%% xyz coordinates (ex: 10 time points)

load('FitData_Batch_log002.mat');
r1 = asdf(1,:)';
r2 = asdf(2,:)';
r3 = asdf(3,:)';

x = zeros(101,3);
y = zeros(101,3);
z = zeros(101,3);

x(:,2) = r2;

x(:,3) = (r1.^2 + r2.^2 - r3.^2)./(2*r2);
y(:,3) = sqrt(r1.^2 - x(:,3).^2);

file_name = "trajectory.pdb";

%% part for file initialization
del = fopen(file_name,"w");
fclose(del);

%% PDB file generation
fileID = fopen(file_name,"a+");

for j=1:100
    fprintf(fileID, "MODEL%9d\n", j);
    for i=1:3
        fprintf(fileID, "%-6s%5d  %-3s%4s %c%4d    %8.3f%8.3f%8.3f%6.2f\n","ATOM", i, "Au", " ", 'A', 0,x(j,i), y(j,i),z(j,i) ,1.0);
    end
    fprintf(fileID,"TER\nENDMDL\n");
end
fclose(fileID);

