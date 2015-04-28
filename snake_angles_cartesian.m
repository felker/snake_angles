%snake_angles_cartesian.m
%Mini-app that solves vacuum RTE in Cartesian coordinates and transforms
%solution to snake coords
%------------------------ PARAMETERS ------------------ %
clear all;
%artificial division:
lx = 7*pi./(2*0.1); %k=0.1 here
ly = lx; 
nx = 40;
ny = nx;

c = 1.0;
dx = lx/nx;
dy = ly/ny;
dt = 0.1;
nt = 80;

N=6;
normalization_tol = 1e-6;
%------------------------ BOUNDARY CONDITIONS  ------------------ %
%%% van Leer interpolation requires two bc 
num_ghost = 2; 
is = 1+num_ghost; ie = nx+num_ghost;
js= 1+num_ghost; je= ny+num_ghost;
%number of physical cells
nx_r= nx;
ny_r = ny;
nx = nx+2*num_ghost;
ny = ny+2*num_ghost;
%------------------------ BUILD MESH ------------------ %
%Radiation Angular Discretization
[ncells,nxa,mu,mu_b,pw] = uniform_angles2D(N);
v = zeros(nx,ny,2); 

%Radiation Spatial Discretization
%Radiation samples are volume averaged quantites
xx=linspace(0,lx,nx_r)';
yy=linspace(0,ly,ny_r)';

%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity = zeros(nx,ny,nxa(1),nxa(2)-1); 

%-------------------- SNAKE COORDINATE DETAILS  ----------------- %
A = 10.0; %amplitude of snake oscillations in y'
K = 0.1; %wavenumber of oscillations DONT USE K AS A LOOP INDEX
beta = A*K*cos(K*xx);
alpha = sqrt(1 + beta.^2);

%--------------- JIANG14 variables, dimensionless constants --------------%
a_r =  5.670373e-8; %stefan boltzman constant
adiabatic = 5/3; %gamma=adiabatic index
R = 1; %8.311; %ideal gas constant

%characteristic values (arbitrary?)
a_0 = 0.1; %characteristic velocity
T_0 = 1.0; %characteristic temperature
P_0 = a_r;%1.0; %characteristic gas pressure

C = c/a_0; %dimensionlesss speed of light
P = a_r*T_0^4/P_0;  %measure of relative importance of rad and gas pressure

%------------------------ PROBLEM SETUP ------------------ %
%Absorption opacities
rho_a = zeros(nx,ny);
%Scattering opacities
rho_s = zeros(nx,ny);
%Fluid density, temperature
density = ones(nx,ny);
temp = ones(nx,ny); 

v(:,:,1) = 0.0*C; 

%------------------------ PRE-TIMESTEPPING SETUP ------------------ %
%Calculate Radiation CFL numbers
%SNAKE EDIT: THIS NEEDS TO BE CALCULATED AT EACH POINT!

% cfl_mu = C*dt*abs(mu)*[1/dx 1/dy]';
% %set dt to have max cfl of 0.4
% dt = dt*1.0/max(cfl_mu);
% %Recalculate radiation CFL
% cfl_mu = C*dt*abs(mu)*[1/dx 1/dy]';
% assert(min(abs(cfl_mu) <= ones(na,1))); 

%------------------------ VELOCITY TERMS ------------------------------ %
[nv, vvnn, vCsquare, vsquare, absV] = update_velocity_terms(v,mu,C);

%------------------------ OUTPUT VARIABLES------------------------------ %
output_interval = 20; 
num_output = 8; %number of data to output
num_pts = nt/output_interval; 
time_out = dt*linspace(0,nt+output_interval,num_pts+1); %extra pt for final step
y_out = zeros(num_pts+1,num_output);

%Explicit-Implicit operator splitting scheme
%-------------------------------------------
for i=0:nt
    time = dt*i; %#ok<NOPTS> 
    %Substep 0: Time series output, boundary condition, moments update
    %Update moments

    if ~mod(i,output_interval)
    end

    %Inject a single ray from the  mu=(1/3,1/3) from center
    %intensity(nx/2,ny/2,1,2) = 1.0; 
    
    %Boundary conditions
    for j=1:num_ghost
        %implement fixed dirichlet absorbing boundary conditions
        intensity(:,js-j,:) = 0.0;
        intensity(:,je+j,:) = 0.0;
        intensity(is-j,:,:) = 0.0;
        intensity(ie+j,:,:) = 0.0;
        %SNAKE periodic boundary conditions TBD
    end
    
    %Inject a ray in the -x +y direction from x_max, y_min (
    for j=1:num_ghost
        intensity(ie+j,js,3,2) = 1.0;
    end
    
    %Substep #1: Explicitly advance transport term
    net_flux = zeros(nx,ny,nxa(1),nxa(2)-1);
    %x-flux
    for j=1:nxa(1) %do all nx, ny at once
        for l=1:nxa(2)-1 %snake edit--- change this discrepancy!
        i_flux = upwind_interpolate2D_snake(mu(j,l,1)*(intensity(:,:,j,l)),dt,dx,ones(nx,ny)*C*sign(mu(j,l,1)),is,ie+1,js,je+1,1);
        net_flux(is:ie,js:je,j,l) = dt*C/dx*(i_flux(is+1:ie+1,js:je) - i_flux(is:ie,js:je));
        end
    end %end of ray trace loop, x-direction
    
    %y-flux
    for j=1:nxa(1) %do all nx, ny at once
        for l=1:nxa(2)-1
        i_flux = upwind_interpolate2D_snake(mu(j,l,2)*(intensity(:,:,j,l)),dt,dy,ones(nx,ny)*C*sign(mu(j,l,2)),is,ie+1,js,je+1,2);
        net_flux(is:ie,js:je,j,l) = net_flux(is:ie,js:je,j,l)+ dt*C/dy*(i_flux(is:ie,js+1:je+1) - i_flux(is:ie,js:je));
        end
    end    %end of ray trace loop, y-direction
    
    intensity = intensity - net_flux; 

    %------------------------ NON-TIME SERIES OUTPUT ------------------ %
    if ~mod(i,output_interval)
        time
       % static_output(); 
        for j=1:nxa(1)
            figure(3);
            %Ray intensity plots
            l=2; %select phi bin
            hi = subplot(2,3,j); 
            h = pcolor(xx,yy,intensity(num_ghost+1:nx_r+2,num_ghost+1:ny_r+2,j,l)');
            %turn off grid lines
            %set(h, 'EdgeColor', 'none');
            x_label = sprintf('mu =(%f, %f)',mu(j,l,1),mu(j,l,2));
            time_title = sprintf('t = %f (s)',time);
            xlabel(x_label);
            title(time_title);
            colorbar
        end
        %Mean intensity plot
         %figure(2);
         %pcolor(xx,yy,rad_energy(num_ghost+1:nx_r+2,num_ghost+1:ny_r+2)')
         pause(0.1)
    end
end

%TRANSFORM INTENSITY SOLUTION to snake angular bins (spatial cells should be the
%same)
    
%COPIED AND CUT DOWN FROM SNAKE ANGLES:
%Renormalize the angular parameterization for snake as a function of
%spatial position
mu_s = zeros(nx_r,ny_r,nxa(1),nxa(2)-1,3);
mu_b_s = zeros(nx_r,ny_r,nxa(1),nxa(2),3);
for i=1:nx_r
    beta_temp = beta(i);
    for j=1:ny_r
        for k=1:nxa(1)
            for l=1:nxa(2)-1
%                 normalization_s = 1./sqrt(1+(beta_temp*mu(k,l,1)).^2 - 2*beta_temp*mu(k,l,1)*mu(k,l,2));
%                 normalization_b = 1./sqrt(1+(beta_temp*mu_b(k,l,1)).^2 - 2*beta_temp*mu_b(k,l,1)*mu_b(k,l,2));
 normalization_s = 1./sqrt((1+beta_temp.^2)*mu(k,l,1).^2 - 2*beta_temp*mu(k,l,1)*mu(k,l,2) + mu(k,l,2).^2 + mu(k,l,3).^2);
 normalization_b = 1./sqrt((1+beta_temp.^2)*mu_b(k,l,1).^2 - 2*beta_temp*mu_b(k,l,1)*mu_b(k,l,2) + mu_b(k,l,2).^2 + mu_b(k,l,3).^2);

 mu_s(i,j,k,l,:) = mu(k,l,:).*normalization_s; 
                mu_b_s(i,j,k,l,:) = mu_b(k,l,:).*normalization_b; 
                %check that each ray has unit spacelike norm in snake coords
                assert(abs(snake_norm(squeeze(mu_s(i,j,k,l,:)),sqrt(1+beta_temp^2),beta_temp) - 1.0) < normalization_tol);
                assert(abs(snake_norm(squeeze(mu_b_s(i,j,k,l,:)),sqrt(1+beta_temp^2),beta_temp)-1.0) < normalization_tol);
            end
            l=nxa(2);
%            normalization_b = 1./sqrt(1+(beta_temp*mu_b(k,l,1)).^2 - 2*beta_temp*mu_b(k,l,1)*mu_b(k,l,2));
 normalization_b = 1./sqrt((1+beta_temp.^2)*mu_b(k,l,1).^2 - 2*beta_temp*mu_b(k,l,1)*mu_b(k,l,2) + mu_b(k,l,2).^2 + mu_b(k,l,3).^2);
    mu_b_s(i,j,k,l,:) = mu_b(k,l,:).*normalization_b; 
            assert(snake_norm(squeeze(mu_b_s(i,j,k,l,:)),sqrt(1+beta_temp^2),beta_temp)==1)
        end
    end
end
%PLOT VARIATION OF ANGLES ALONG X
l=2; %mu(:,2,3)= 0.3827
%at a particular \phi' 
figure(10);
for n=nx_r:-1:1
    p= 1;
    beta_temp = beta(n);
    for k=1:nxa(1)
        %Transform snake rays to cartesian basis, DONT RENORMALIZE
        mu_unprimed = zeros(3,1);
        mu_unprimed(1) = mu_s(n,p,k,l,1);
        mu_unprimed(2) = (mu_s(n,p,k,l,2) - A*K*cos(K*xx(n))*mu_s(n,p,k,l,1));
        mu_unprimed(3) = mu_s(n,p,k,l,3);
        %check normalization of vector in cartesian basis
        assert(abs(mu_unprimed'*mu_unprimed - 1.0) < normalization_tol);
        theta_s = atan2(mu_unprimed(2),mu_unprimed(1)); 
        quiver(0,0,cos(theta_s),sin(theta_s),0,'-r');
        hold on; 
        axis([-1 1 -1 1]);
        mu_b_unprimed = zeros(3,1);
        mu_b_unprimed(1) = mu_b_s(n,p,k,l,1);
        mu_b_unprimed(2) = (mu_b_s(n,p,k,l,2) - A*K*cos(K*xx(n))*mu_b_s(n,p,k,l,1));
        mu_b_unprimed(3) = mu_b_s(n,p,k,l,3);
        assert(abs(mu_b_unprimed'*mu_b_unprimed - 1.0) < normalization_tol);
        theta_b = atan2(mu_b_unprimed(2),mu_b_unprimed(1)); 
        quiver(0,0,cos(theta_b),sin(theta_b),0,'-k');        
        titlestr = ['x =',num2str(xx(n)),',\mbox{   }','$$\frac{2xk}{\pi} =$$',num2str(xx(n)*2*K/pi)];
        title(titlestr,'Interpreter','latex');
        %quiver3(0,0,0,mu_s(n,p,k,l,1),mu_s(n,p,k,l,2),mu_s(n,p,k,l,3))
    end
    hold off; 
    pause(); 
end

count = 0; 
intensity_snake = zeros(nx_r,ny_r,nxa(1),nxa(2)-1); 
for n=1:nx_r %should use overlap of solid angle bins
    for p=1:ny_r
        for j=1:nxa(1) % even though we only have one nonzero theta bin in current 
            %cartesian problem
            for l=2:2
                %Transform from Cartesian basis to Snake Basis
                mu_prime = zeros(3,1); 
                mu_prime(1) = mu(j,l,1);
                mu_prime(2) = (mu(j,l,2) + A*K*cos(K*xx(n))*mu(j,l,1));
                mu_prime(3) = mu(j,l,3);
                %check that transformed vector has norm 1 in snake coords
                assert(abs(snake_norm(mu_prime,sqrt(1+beta(n)^2),beta(n)) - 1.0) < normalization_tol);

%           THE ISSUE LIES WITH THE SNAKE RENORMALIZATION MESSING WITH THE
%           MU_Z COMPONENTS

                %debug at one point in space
%                 if (n==nx_r-2 && p== 4)
%                     j
%                     intensity(n,p,j,l)
%                     mu(j,l,:)
%                     mu_s(n,p,j,l,:)
%                     mu_prime
%                 end
                %find new mu bin
                for k=1:nxa(1)
                    if k==nxa(1) %theta bins are ordered and periodic
                        neighbor = 1; 
                    else
                        neighbor = k+1;
                    end
                      %sort by theta bin! use 4 quadrant arctangent
                      transformed_theta_prime = atan2(mu_prime(2),mu_prime(1)) +pi;
                      lower_b_theta = atan2(mu_b_s(n,p,k,l,2),mu_b_s(n,p,k,l,1)) + pi; 
                      upper_b_theta = atan2(mu_b_s(n,p,neighbor,l,2),mu_b_s(n,p,neighbor,l,1)) +pi;
                      if (upper_b_theta == 0.0) %this happens a lot
                          upper_b_theta = 2*pi; 
                      end
                      
                      if (transformed_theta_prime >= lower_b_theta) && (transformed_theta_prime ...
                              <= upper_b_theta)
                        intensity_snake(n,p,k,l) = intensity_snake(n,p,k,l) +...
                            intensity(n+num_ghost,p+num_ghost,j,l);
                        count = count+1;                         
                      end
                end                         
            end
        end
    end
end
count %since current algorithm maps 1-1 from cartesian thetas to
%snake thetas at each (nx)(ny) cell
%I have checked that sum(sum(sum(sum(intensity_snake)))) = sum(sum(sum(sum(intensity))))
%after taking out BCs
 for j=1:nxa(1)
            figure(4);
            %Ray intensity plots
            l=2; %select phi bin
            hi = subplot(2,3,j); 
            h = pcolor(xx*ones(1,nx_r),ones(ny_r,1)*yy'+A*sin(xx*ones(1,nx_r).*K),intensity_snake(:,:,j,l));
            %turn off grid lines
            %set(h, 'EdgeColor', 'none');
            x_label = sprintf('mu =(%f, %f)',mu(j,l,1),mu(j,l,2));
            time_title = sprintf('t = %f (s)',time);
            xlabel(x_label);
            title(time_title);
            colorbar
 end
 