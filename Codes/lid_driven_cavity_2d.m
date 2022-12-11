tStart = cputime;
% Constant parameters
nx = 128;
ny = 128;
nu = 0.001; % dynamic viscosity
lx = 1;
ly = 1;
dx = lx/nx;
dy = ly/ny;
%% Boundary conditions
% All normal velocity components are zero - no inlets whatsoever
uT = 1.0;
uB = 0.0;
vL = 0.0;
vR = 0.0;
% Reynolds no.
Re = uT*lx/nu;
% Choose dt based on linear advection-diffusion constraint
C = 1;
dt_max = min([dx/uT,2^(2/3)*C^(1/3)*(dx/uT)^(4/3) , 0.25*dx*dx/nu]);
%% CFL condition
dt = dt_max;
CFL_no = uT*dt/dx;
% dt = 1e-6;
%% Initialize variables - we will consider a staggered grid to the minus side
% All variables will be of size (nx + 2, ny + 2)
u = zeros(ny+2,nx+2);
v = zeros(ny+2,nx+2);
u_old = u;
v_old = v;
% usol = [];
% usol = [usol,u];
% vsol = [];
% vsol = [vsol,v];
change_in_u = zeros(ny+2,nx+2);
change_in_v = zeros(ny+2,nx+2);
l2_norm_u = [];
l2_norm_v = [];
err_u = 1;
err_v = 1;
ut = zeros(ny+2,nx+2);
vt = zeros(ny+2,nx+2);
dummy_convection_x = zeros(ny+2,nx+2);
dummy_diffusion_x = zeros(ny+2,nx+2);
dummy_convection_y = zeros(ny+2,nx+2);
dummy_diffusion_y = zeros(ny+2,nx+2);
p = zeros(ny+2,nx+2); 
% psol = [];
% psol = [psol,p];
convection = zeros(1,1);
diffusion = zeros(1,1);
%% Build the pressure coefficients
Ae = 1/dx^2*ones(ny+2,nx+2);
Aw = 1/dx^2*ones(ny+2,nx+2);
An = 1/dy^2*ones(ny+2,nx+2);
As = 1/dy^2*ones(ny+2,nx+2);
% set the left wall coefficients
Aw(2:end-1,2) = 0.0;
% set the right wall coefficients
Ae(2:end-1,end-1) = 0.0;
% set the top wall coefficients
An(2,2:end-1) = 0.0;
% set the bottom wall coefficients
As(end-1,2:end-1) = 0.0;
Ap = -(Aw + Ae + An + As);
%% Time integration
n = 1;
while err_u > 1e-8 && err_v > 1e-8
% while n < 100
    % Set bcs on u velocity
    % left wall
    u(2:end-1,2) = 0.0;
    % right wall
    u(2:end-1,end) = 0.0;
    % top wall
    u(1,2:end) = 2.0*uT - u(2,2:end);
    % bottom wall
    u(end,2:end) = 2.0*uB - u(end-1,2:end);
    % Set bcs on v velocity
    % top wall
    v(1,2:end-1) = 0.0;
    % bottom wall
    v(end-1,2:end-1) = 0.0;
    % left wall
    v(1:end-1,1) = 2.0*vL - v(1:end-1,2);
    % right wall
    v(1:end-1,end) = 2.0*vR - v(1:end-1,end-1);
%% Assemble x-mom first for u_tilde - do interior points only
    for i = 3:nx+1
        for j = 2:ny+1
            ue = 0.5*(u(j,i)+u(j,i+1));
            uw = 0.5*(u(j,i)+u(j,i-1));
            un = 0.5*(u(j,i)+u(j+1,i));
            us = 0.5*(u(j,i)+u(j-1,i));
            vn = 0.5*(v(j+1,i-1)+v(j+1,i));
            vs = 0.5*(v(j,i)+v(j,i-1));
            convection(j,i) = -(ue*ue - uw*uw)/dx - (un*vn - us*vs)/dy;
            diffusion(j,i) = (1/Re)*((u(j,i+1) - 2*u(j,i) + u(j,i-1))/dx^2 + ...
                        (u(j+1,i) - 2*u(j,i) + u(j-1,i))/dy^2) ;
            ut(j,i) = u(j,i) + dt/2*(3*(convection(j,i) + diffusion(j,i)) - (dummy_convection_x(j,i) + dummy_diffusion_x(j,i))); % 2nd order Adams Bashforth
%             ut(j,i) = u(j,i) + dt*(convection(j,i) + diffusion(j,i));
            dummy_convection_x(j,i) = convection(j,i);
            dummy_diffusion_x(j,i) = diffusion(j,i);
        end
    end
    %% Assemble y-mom first for v_tilde - do interior points only
    for i = 2:nx+1
        for j = 2:ny
            ue = 0.5*(u(j,i+1) + u(j-1,i+1));
            uw = 0.5*(u(j-1,i) + u(j,i));
            ve = 0.5*(v(j,i) + v(j,i+1));
            vw = 0.5*(v(j,i) + v(j,i-1));
            vn = 0.5*(v(j,i) + v(j+1,i));
            vs = 0.5*(v(j,i) + v(j-1,i));
            convection(j,i) = -(ve*ue - vw*uw)/dx - (vn*vn - vs*vs)/dy;
            diffusion(j,i) = (1/Re)*((v(j,i+1) - 2*v(j,i) + v(j,i-1))/dx^2 + ...
                        (v(j+1,i) - 2*v(j,i) + v(j-1,i))/dy^2) ;
            vt(j,i) = v(j,i) + dt/2*(3*(convection(j,i) + diffusion(j,i)) - (dummy_convection_y(j,i) + dummy_diffusion_y(j,i))); % 2nd order Adams Bashforth
%             vt(j,i) = v(j,i) + dt*(convection + diffusion);
            dummy_convection_y(j,i) = convection(j,i);
            dummy_diffusion_y(j,i) = diffusion(j,i);
        end
    end
    %% Compute pressure rhs: prhs = 1/dt * div(ut)
    divut = zeros(ny+2,nx+2);
    divut(2:end-1,2:end-1) = (ut(2:end-1,3:end) - ut(2:end-1,2:end-1))/dx + (vt(2:end-1,2:end-1) - vt(1:end-2,2:end-1))/dy;
    prhs = divut/dt;
    %% Solve pressure poisson equation
    p = sor_solver(p,Ae,Aw,An,As,Ap,prhs,nx,ny);
    % Final time advance - do interior only
    % u = ut - dt*dpdx
    u(2:end-1,3:end-1) = ut(2:end-1,3:end-1) - dt*(p(2:end-1,3:end-1)-p(2:end-1,2:end-2))/dx; % (n+1)th time step
    v(2:end-2,2:end-1) = vt(2:end-2,2:end-1) - dt*(p(3:end-1,2:end-1)-p(2:end-2,2:end-1))/dy; % (n+1)th time step
%     usol = [usol,u];
%     vsol = [vsol,u];
    change_in_u = u - u_old;
    change_in_v = v - v_old;
%     disp(size(change_in_u(:,:,n)));
    l2_norm_u = [l2_norm_u,sqrt(sum(change_in_u.^2,"all")/(nx*ny))];
    l2_norm_v = [l2_norm_v,sqrt(sum(change_in_v.^2,"all")/(nx*ny))];
    err_u = l2_norm_u(end);
    err_v = l2_norm_v(end);
    fprintf('\t%d \t%f \t%f\n',n,l2_norm_u,l2_norm_v);
    n = n+1;
    u_old = u;
    v_old = v;
end
x_coord = dt*linspace(2,n,n-1);
figure(1)
xlim([0,80])
% ylim([1e-8,1])
% subplot(2,1,1)
title('Re = 100 and Re = 1000 on a 128*128 finite volume grid')
semilogy(x_coord,l2_norm_u,'^');
xlabel('non-dimensional time units')
ylabel('average l2 norm of u vel')
hold on
legend('Re = 100,\delta t = 1.5e-3','Re = 1000, \delta t = 2.46e-3')
% subplot(2,1,2);
% semilogy(x_coord,l2_norm_v);
% xlabel('non-dimensional time units')
% ylabel('average l2 norm of v vel')
% figure(2)
% image(divut);
% colorbar
tEnd = cputime - tStart;
y_coord = [0,0.0547,0.0625,0.0703,0.1016,0.1719,0.2813,0.4531,0.5,0.6172,0.7344,0.8516,0.9531,0.9609,0.9688,0.9766,1]';
x_coord = [0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1]';
u_ghia_Re100 = [0,-0.03717,-0.04192,-0.04775,-0.06434,-0.10150,-0.15662,-0.21090,-0.20581,-0.13641,0.00332,0.23151,0.68717,0.73722,0.78871,0.84123,1]';
v_ghia_Re100 = [0,0.09233,0.10091,0.10890,0.12317,0.16077,0.17507,0.17527,0.05454,-0.24533,-0.22445,-0.16914,-0.10313,-0.08864,-0.07391,-0.05906,0]';
u_ghia_Re1000 = [0,-0.18109,-0.20196,-0.22220,-0.29730,-0.38289,-0.27805,-0.10648,-0.06080,0.05702,0.18719,0.33304,0.46604,0.51117,0.57492,0.65928,1]';
v_ghia_Re1000 = [0,0.27485,0.29012,0.30353,0.32627,0.37095,0.33075,0.32235,0.02526,-0.31966,-0.42665,-0.51550,-0.39188,-0.33714,-0.27669,-0.21388,0]';
figure(2)
title('u velocity along vertical line through geometric centre of cavity')
plot(y_coord,u_ghia_Re1000)
hold on
% % plot(linspace(0,1,66),u(end:-1:1,34),'--')
% % hold on
plot(linspace(0,1,130),u(end:-1:1,66),'o')
% % hold on
% plot(linspace(0,1,34),u(end:-1:1,18),'d')
figure(3)
plot(x_coord,v_ghia_Re1000)
hold on
% plot(linspace(0,1,66),v(33,end:-1:1),'--')
% hold on
plot(linspace(0,1,130),v(65,end:-1:1),'o')
% % hold on
% plot(linspace(0,1,34),v(17,end:-1:1),'d')
% %% n * n finite volumes
% y_grid_pt = round(1 + (1-y_coord)*(nx+1)); % at x = 0.5
% x_grid_pt = round(1 + x_coord*(nx+1)); % at y = 0.5
% 
% % Book-keeping of u and v velocities at above grid points
% u_midway = zeros(1,1);
% v_midway = zeros(1,1);
% for i = 1:length(y_grid_pt) 
%     u_midway(i,1) = u(y_grid_pt(i),nx/2+2);
% end
% for i = 1:length(x_grid_pt)
%     v_midway(i,1) = v(nx/2+1,x_grid_pt(i));
% end
% figure(2)
% plot(y_coord,u_ghia_Re100,y_coord,u_midway(:,1),'o',y_coord,u_midway(:,2),'s',y_coord,u_midway(:,3),'^')
% xlabel('y coordinate at x = 0.5')
% ylabel('u velocity')
% legend('Ghia paper Table I','32*32 grid,\delta t = 3.125e-3','64*64 grid,\delta t = 5e-4','128*128 grid,\delta t = 1.25e-4');
% figure(3)
% plot(x_coord,v_ghia_Re100,x_coord,v_midway(:,1),'o',x_coord,v_midway(:,2),'s',x_coord,v_midway(:,3),'^')
% legend('Ghia paper Table II','32*32 grid,\delta t = 3.125e-3','64*64 grid,\delta t = 5e-4','128*128 grid,\delta t = 1.25e-4');
%% SOR solver function
function phi = sor_solver(phi_0,Ae,Aw,An,As,Ap,rhs_vec,nx,ny)
   w = 2/(1+sin(pi/(nx+1)));
   phi = phi_0;
   r = zeros(nx,ny);
   residual = 1;
%    error = 1;
   tol = 1e-5;
%    residue_history = [];
%    error_history = [];
%    phi = zeros(nx+2,ny+2);
   while residual > tol
       for i = 2:nx+1
           for j = 2:ny+1
               phi(j,i) = w*((rhs_vec(j,i) - Ae(j,i)*phi(j,i+1) - Aw(j,i)*phi(j,i-1) - ...
                   An(j,i)*phi(j+1,i) - As(j,i)*phi(j-1,i))/Ap(j,i)) + (1-w)*phi(j,i);
           end
       end
       residual = 0;
       for i = 2:nx+1
           for j = 2:ny+1
               r(j,i) = rhs_vec(j,i) - Ae(j,i)*phi(j,i+1) - Aw(j,i)*phi(j,i-1) - ...
                   An(j,i)*phi(j+1,i) - As(j,i)*phi(j-1,i) - Ap(j,i)*phi(j,i);
               residual = residual + r(j,i)*r(j,i);
           end
       end
       residual = sqrt(residual/(nx*ny));
%        error = sqrt(sum((phi - phi_0).^2,"all")/(nx*ny));
%        residue_history = [residue_history,residual];
%        error_history = [error_history,sqrt(sum(error,'all')/M^2)];
%        disp(residual);
   end
end
