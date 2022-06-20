function [M] = grid_laplace(L,Nx,BC,D)
%GRID-LAPLACE 均匀正方网格拉普拉斯算子
x = linspace(0,L,Nx);
dx = x(2)-x(1);


%initialize Matrices
 N = Nx^2;
 M = zeros(N,N);
 
 
 %内点
 for i =2:Nx-1
     for j = 2:Nx-1
         n=i+(j-1)*Nx;
         M(n,n) = -4;
         M(n,n-1) = 1;
         M(n,n+1) = 1;
         M(n,n-Nx) = 1;
         M(n,n+Nx) = 1;
     end
 end
  
 


 %down and up 内点
 n = (Nx-1)*Nx;
 for j = 2:Nx-1
     M(j,j) = -2;
     M(j,j-1)=1;
     M(j,j+1)=1;
     M(n+j,n+j)=-2;
     M(n+j,n+j-1)=1;
     M(n+j,n+j+1)=1;
     
 end
   
 %left and Right 内点
 for i = 2:Nx-1
     M(Nx*(i-1)+1,Nx*(i-1)+1)=-2;
     M(Nx*(i-1)+1,Nx*(i-2)+1)=1;
     M(Nx*(i-1)+1,Nx*i+1)=1;
     M(Nx*i,Nx*i)=-2;
     M(Nx*i,Nx*(i-1))=1;
     M(Nx*i,Nx*(i+1))=1;
 end


 
 %down and up 内点
 for j = 2:Nx-1
    [f1,f2] = Force_ring_xy(x(j),0);
     M(j,j) = -2*f2*dx/D-2;
     M(j,j+Nx) = 2;
    [f1,f2] = Force_ring_xy(x(j),L);
     M(n+j,n+j)=2*f2*dx/D+2;
     M(n+j,n+j-Nx)=-2;
 end
 
  %left and Right 内点
 for i = 2:Nx-1
     [f1,f2] = Force_ring_xy(0,x(i));
     M(Nx*(i-1)+1,Nx*(i-1)+1)=-2*f1*dx/D-2;
     M(Nx*(i-1)+1,Nx*(i-1)+2)=2;
     [f1,f2] = Force_ring_xy(L,x(i));
     M(Nx*i,Nx*i)=2*f1*dx/D+2;
     M(Nx*i,Nx*i-1)=-2;  
 end
 
 %corner
 %左上
 [f1,f2] = Force_ring_xy(0,L);
 n=Nx*(Nx-1);
 M(n+1,n+1) = 2*(f2-f1)*D*dx -4;
 M(n+1,n-Nx+1) = 2;
 M(n+1,n+2)=2;
 
 %右上
 [f1,f2] = Force_ring_xy(L,L);
 n = Nx^2;
 M(n,n) = 2*dx*(f1+f2)/D-4;
 M(n,n-1)=2;
 M(n,n-Nx)=2;
 
 %左下
 [f1,f2] = Force_ring_xy(0,0);
 M(1,1)=-4-2*dx*(f1+f2)/D;
 M(1,2)=2;
 M(1,1+Nx)=2;
 
 %右下
 [f1,f2] = Force_ring_xy(L,0);
 M(Nx,Nx) = 2*(f1-f2)*dx/D-4;
 M(Nx,Nx-1)=2;
 M(Nx,2*Nx)=2;
 
 

 

 if BC == 0 %Dirichlet边界 J=0
 %up and down 内点
 n = (Nx-1)*Nx;
 for j = 2:Nx-1

     
 end
   
 %left and Right 内点
 for i = 2:Nx-1
  
 end
 
 

    
    
    


 
 
end


