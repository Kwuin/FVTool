Lx = 4;
Ly = 4;
Nx = 100;
Ny = 100;
m = createMesh2D(Nx,Ny,Lx,Ly);


%边界
D_val = 0.005; %diffusion coeeficient
BC = createBC(m); % all Neumann boundary condition structure
a = BC.left.a;
BC.left.a = -D_val*a;
BC.right.a = -D_val*a;
BC.top.a = -D_val*a;
BC.bottom.a = -D_val*a;

[xb,yb] = Force_diff_xy(0,m.cellcenters.y); %diff/ring
BC.left.b = transpose(xb);
[xb,yb] = Force_diff_xy(Lx,m.cellcenters.y);
BC.right.b = transpose(xb);

[xb,yb] = Force_diff_xy(m.cellcenters.x,0);
BC.bottom.b = transpose(yb);
[xb,yb] = Force_diff_xy(m.cellcenters.x,Ly);
BC.top.b = transpose(yb);

c = 0.01*ones(1,Nx);
BC.left.c = -c;
BC.right.c = c;
BC.top.c = c;
BC.bottom.c = -c;

%diffussion term
D = createCellVariable(m, D_val); % assign dif. coef. to all the cells
Dave = harmonicMean(D); % convert a cell variable to face variable

%convection term
u_face = createFaceVariable(m, 1); % assign velocity value to cell faces
[x,y] = ndgrid(m.facecenters.x, m.cellcenters.y);
[xv,yv] = Force_diff_xy(x,y);
u_face.xvalue = xv;

[x,y] = ndgrid(m.cellcenters.x, m.facecenters.y)
[xv,yv] = Force_diff_xy(x,y);
u_face.yvalue = yv;

%solve
Mconv =  convectionTerm(u_face); % convection term, central, second order
Mdiff = diffusionTerm(Dave); % diffusion term
[Mbc, RHS] = boundaryCondition(BC); % boundary condition discretization
M = Mconv-Mdiff+Mbc; % matrix of coefficient for central scheme
c = solvePDE(m, M, RHS); % solve for the central scheme
c.value = c.value-min(min(c.value));
c.value = log(c.value);
%c.value = 
visualizeCells2D(c);














