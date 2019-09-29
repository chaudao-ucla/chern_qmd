var1 = importdata('var1.dat');
var2 = importdata('var2.dat');
var3 = importdata('var3.dat');
var4 = importdata('var4.dat');
var5 = importdata('var5.dat');
var6 = importdata('var6.dat');
var7 = importdata('var7.dat');
var8 = importdata('r1.dat');
myvar = importdata('r1.dat');
numvar = importdata('numvar.dat');

plot(table2array(r1(:,1)),table2array(r1(:,2)));

%pos1 = importdata('pos1.dat');
%scatter(pos1(:,1),pos1(:,2));

%plot(var8(:,1),var8(:,2));

%{
plot(var2(:,1),var2(:,2)*2.0,var2(:,1),var2(:,2),var3(:,1),var3(:,2),var4(:,1),var4(:,2),var5(:,1),var5(:,2),var6(:,1),var6(:,2),var7(:,1),var7(:,2));
title('Energy Variance for dt = 0.001, N = 30, Num\_steps = 2000');
ylabel('Variance')
xlabel('Density')
legend('r=1','r=2','r=3','r=4','r=5','r=6','r=7');
%}
%{
plot(numvar(:,1), numvar(:,2))
title('Energy variance for dt = 0.001, r=5, Num\_steps = 5000');
ylabel('Variance');

xlabel('N');
%}

%{
plot(var3(:,1),var3(:,2),var4(:,1),var4(:,2),var5(:,1),var5(:,2),var6(:,1),var6(:,2),var7(:,1),var7(:,2));
title('Energy Variance for dt = 0.001, N = 30, Num\_steps = 2000');
ylabel('Variance')
xlabel('Density')
legend('r=3','r=4','r=5','r=6','r=7');
%}

%{
plot(var2(:,1)/1.0,var2(:,2)*2.0,var2(:,1)/2.0,var2(:,2),var3(:,1)/3.0,var3(:,2),var4(:,1)/4.0,var4(:,2),var5(:,1)/5.0,var5(:,2),var6(:,1)/6.0,var6(:,2),var7(:,1)/7.0,var7(:,2));
title('Energy Variance for dt = 0.001, N = 30, Num\_steps = 2000');
ylabel('Variance')
xlabel('Density/Cutoff Radius')
legend('r=1','r=2','r=3','r=4','r=5','r=6','r=7');
%}

