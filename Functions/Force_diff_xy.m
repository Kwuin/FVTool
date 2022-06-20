function [f1,f2] = Force_diff_xy(x1,x2)

%分化图景的动力对坐标矩阵的函数

a =1;
b=1;
S=0.5;
n=4;
k=1;
Sn = S^n;
f1 = a*x1.^n./(Sn+x1.^n) + b*S^n./(Sn+x2.^n)-k*x1;
f2 = a*x2.^n./(Sn+x2.^n) + b*S^n./(Sn+x1.^n)-k*x2;
end