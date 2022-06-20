function [f1,f2] = Force_ring_xy(X,Y)
%环流动力-坐标juzhen函数  X Y分别为网格的横纵坐标矩阵
a = 0.1;
b = 0.1;
c = 100;
eps = 0.1;
tau = 5;
f1 = (eps^2+X.^2)./((1+X.^2).*(1+Y))-a*X;
f2 = (b-Y./(1+c*X.^2))/tau;
end

