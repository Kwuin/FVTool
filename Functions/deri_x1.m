function [M] = deri_x1(L,Nx,D,bc,f)
%Force = ring
x = linspace(0,L,Nx);
dx = x(2)-x(1);

%DERI_X x方向一阶导数
N = Nx^2;
M = zeros(N,N);

%内点
for i= 2:Nx-1
    for j = 2:Nx-1
        n = Nx*(i-1);
        M(n+j,n+j+1)=1/2;
        M(n+j,n+j-1)=-1/2;
    end
end
M = M/dx;


%left and right
for i = 1:Nx
    [f1,f2]=Force_ring_xy(0,(i-1)*dx);
    M((i-1)*Nx+1,(i-1)*Nx+1) = f1/D;
    [f1,f2]=Force_ring_xy(L,(i-1)*dx);
    M(i*Nx, i*Nx) = f1/D;
end

if bc == 1
   n=Nx*(Nx-1);
   for j=2:Nx-1
       [f1,f2]=Force_ring_xy((j-1)*dx,0);
       M(j,j)=f1/D;
       [f1,f2]=Force_ring_xy((j-1)*dx,L);
       M(n+j,n+j)=f1/D;
   end
else
    for j=2:Nx-1
       M(j,j-1)=-1/(2*dx);
       M(j,j+1)=1/(2*dx);
       M(n+j,n+j-1)=-1/(2*dx);
       M(n+j,n+j+1)=1/(2*dx);
   end
    
    
    
end
end
