%snake_angles.m
%Mini-app that solves vacuum RTE. we inject an initial propagation of rays
%and then transform these rays to snake coordinates to see how the angular
%mesh evolves. 

%this only tests advection step

%------------------------ PARAMETERS ------------------ %
clear all;
%close all;
% nx = 100;
% ny = 100;
% lx = 100.0;
% ly = 100.0; %corresponding to max y_cartesian, not the actual max y' 
N = 6; 

%artificial division:
lx = 7*pi./(2*0.1); %k=0.1 here
ly = lx; 
nx = 40;
ny = nx;

c = 1.0;
dx = lx/nx;
dy = ly/ny;
dt = 0.1;
nt = 200;
normalization_tol = 1e-6;

%------------------------ BOUNDARY CONDITIONS  ------------------ %
%%% van Leer interpolation requires two bc 
num_ghost = 2; 
is = 1+num_ghost; ie = nx+num_ghost;
js= 1+num_ghost; je= ny+num_ghost;
%number of physical cells
nx_r = nx;
ny_r = ny;
nx = nx+2*num_ghost;
ny = ny+2*num_ghost;
%------------------------ BUILD SPATIAL MESH ------------------ %
%Radiation Spatial Discretization
%Radiation samples are volume averaged quantites
xx=linspace(0,lx,nx_r)';
yy=linspace(0,ly,ny_r)';

%-------------------- SNAKE COORDINATE DETAILS  ----------------- %
A = 10.0; %amplitude of snake oscillations in y'
K = 0.1; %wavenumber of oscillations DONT USE K AS A LOOP INDEX... too late
beta = A*K*cos(K*xx);
alpha = sqrt(1 + beta.^2);

%------------------------ BUILD ANGULAR MESH ------------------ %
%Radiation Angular Discretization
%This produces a local discretization when the metric is in the canonical
%form
[ncells,nxa,mu,mu_b,pw] = uniform_angles2D(N);
%Renormalize the angular parameterization for snake as a function of
%spatial position
% \hat{k}^i' (x^\mu) = CW^i
%nx_r or nx?
mu_s = zeros(nx,ny,nxa(1),nxa(2)-1,3);
mu_b_s = zeros(nx,ny,nxa(1),nxa(2),3);
for i=1:nx
    if i >= is && i<=ie %deal with ghost cells for which beta is undefined
        beta_temp = beta(i-num_ghost); 
    elseif i < is
        beta_temp = beta(1); 
    elseif i > ie
        beta_temp = beta(nx_r); 
       % beta_temp = 1.0; %this might screw up the van leer part....
        %Yes, this will cause problems in nonzero dirichlet bcs and
        %periodic bcs. but if we are setting absorbing boundaries, then it
        %shouldnt matter. 
        
        % for now, extend the mu arrays at the boundaries by copying them
        % to the ghosts
    end
    for j=1:ny
        for k=1:nxa(1)
            for l=1:nxa(2)-1
                normalization_s = 1./sqrt(1+(beta_temp*mu(k,l,1)).^2 - 2*beta_temp*mu(k,l,1)*mu(k,l,2));
                normalization_b = 1./sqrt(1+(beta_temp*mu_b(k,l,1)).^2 - 2*beta_temp*mu_b(k,l,1)*mu_b(k,l,2));
                mu_s(i,j,k,l,:) = mu(k,l,:).*normalization_s; 
                mu_b_s(i,j,k,l,:) = mu_b(k,l,:).*normalization_b; 
                %check that each ray has unit spacelike norm in snake coords
                assert(abs(snake_norm(squeeze(mu_s(i,j,k,l,:)),sqrt(1+beta_temp^2),beta_temp) - 1.0) < normalization_tol);
                assert(abs(snake_norm(squeeze(mu_b_s(i,j,k,l,:)),sqrt(1+beta_temp^2),beta_temp)-1.0) < normalization_tol);
            end
            l=nxa(2);
            normalization_b = 1./sqrt(1+(beta_temp*mu_b(k,l,1)).^2 - 2*beta_temp*mu_b(k,l,1)*mu_b(k,l,2));
            mu_b_s(i,j,k,l,:) = mu_b(k,l,:).*normalization_b; 
            assert(snake_norm(squeeze(mu_b_s(i,j,k,l,:)),sqrt(1+beta_temp^2),beta_temp)==1)
        end
    end
end

%------------------------ INTENSITY AND FLUID VELOCITY ------------------ %
%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity = zeros(nx,ny,nxa(1),nxa(2)-1); 
v = zeros(nx,ny,2); 

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
% CFL CONDITION IS NOT TRIVIAL-- need to do stability analysis

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
    %No! inject the ray from x_max, y_middle to avoid corner effects
    for j=1:num_ghost
        intensity(ie+j,js+ny/2,3,2) = 1.0;
    end
    
    %Substep #1: Explicitly advance transport term
    net_flux = zeros(nx,ny,nxa(1),nxa(2)-1);
    %x-flux
    for j=1:nxa(1) %do all nx, ny at once
        for l=1:nxa(2)-1 %snake edit: with new normalized \hat{k}, how will it work?
            %cannot pull mu out from partial, since it changes with x
            %position
        %i_flux = upwind_interpolate2D_snake(mu(j,l,1)*(intensity(:,:,j,l)),dt,dx,ones(nx,ny)*C*sign(mu(j,l,1)),is,ie+1,js,je+1,1);
        i_flux = upwind_interpolate2D_snake(mu_s(:,:,j,l,1).*(intensity(:,:,j,l)),dt,dx,C.*sign(mu_s(:,:,j,l,1)),is,ie+1,js,je+1,1);
        net_flux(is:ie,js:je,j,l) = dt*C/dx*(i_flux(is+1:ie+1,js:je) - i_flux(is:ie,js:je));
        end
    end %end of ray trace loop, x-direction
    
    %y-flux
    for j=1:nxa(1) %do all nx, ny at once
        for l=1:nxa(2)-1
        %i_flux = upwind_interpolate2D_snake(mu(j,l,2)*(intensity(:,:,j,l)),dt,dy,ones(nx,ny)*C*sign(mu(j,l,2)),is,ie+1,js,je+1,2);
        i_flux = upwind_interpolate2D_snake(mu_s(:,:,j,l,2).*(intensity(:,:,j,l)),dt,dy,C.*sign(mu_s(:,:,j,l,2)),is,ie+1,js,je+1,2);
        net_flux(is:ie,js:je,j,l) = net_flux(is:ie,js:je,j,l)+ dt*C/dy*(i_flux(is:ie,js+1:je+1) - i_flux(is:ie,js:je));
        end
    end    %end of ray trace loop, y-direction
    
    %Substep #1.1: Compute solid angular fluxes
    %the coefficients C^i can depend C^i(x,y,z,k^1,k^2,k^3), so the outer
    %loop should be over solid angle, and the inner loops will be over
    %spatial cells... actually the order shouldnt matter
        
    %Compute C^2
    %for this metric/coords, this can be precomputed. varies from -|C2| to
    %|C2| due to sin varying over entire range on this mesh
     C2 = zeros(nx_r,ny_r,nxa(1),nxa(2)-1);
     for n=1:nx_r 
         for p=1:ny_r
             for j=1:nxa(1)
                 for l=1:nxa(2)-1
                    %spatial dependence: beta(x), \hat{k}^1'(x) from normalization
                    C2(n,p,j,l) = A*K^2*sin(K*xx(n))*(1+2*beta(n).^2)*mu_s(n+num_ghost,p+num_ghost,j,l,1).^2;              
                 end
             end
         end
     end
     %Compute parameterization speed
     dxadk = zeros(nx_r,ny_r,nxa(1),nxa(2)-1); %this shouldnt depend on \hat{k}^i' = mu_s!!
     for n=1:nx_r 
         for p=1:ny_r
             for j=1:nxa(1)
                 for l=1:nxa(2)-1
                     temp = 1+(beta(n)*mu(j,l,1)).^2 - 2*beta(n)*...
                         mu(j,l,1)*mu(j,l,2);
                     dxadk(n,p,j,l) = beta(n)*mu(j,l,1).^2./(temp)^(3/2)*sqrt(1- ...
                        mu(j,l,1).^2./(temp));
                 end
             end
         end
     end
     %Manually hardcode the donor cell method since we only have vtheta
     angular_flux = zeros(nx_r,ny_r,nxa(1),nxa(2)-1);
     i_flux = zeros(nx_r,ny_r,nxa(1),nxa(2)-1,2);
     dtheta = 2*pi./(nxa(1)); 
     dphi = pi./(nxa(2));
     for n=1:nx_r
        for p=1:ny_r
            for j=1:nxa(1) %edit out for easy periodic circular shifts
                for l=1:nxa(2)-1
                    v1 = C2(n,p,j,l).*dxadk(n,p,j,l);
                    %theta flux
                    if v1 > 0 %counterclockwise speed
                        if j-1 == 0
                          i_flux(n,p,j,l,1) = intensity(n+num_ghost,p+num_ghost,nxa(1),l,1); %flux entering "left" boundary                          
                        else
                          i_flux(n,p,j,l,1) = intensity(n+num_ghost,p+num_ghost,j-1,l,1); %flux entering "left" boundary                          
                        end
                    elseif v1 < 0 %clockwise 
                        i_flux(n,p,j,l,1) = intensity(n+num_ghost,p+num_ghost,j,l,1);  %flux leaving left boundary          
                    end
                end
                %phi flux
                %handle both poles
            end %end of angular loops
            for j=1:nxa(1)
                for l=1:nxa(2)-1
                    v1 = C2(n,p,j,l).*dxadk(n,p,j,l);
                    if j+1 == nxa(1)+1
                        angular_flux(n,p,j,l) = dt*v1/dtheta.*(i_flux(n,p,1,l,1) - i_flux(n,p,j,l,1));                   
                    else
                        angular_flux(n,p,j,l) = dt*v1/dtheta.*(i_flux(n,p,j+1,l,1) - i_flux(n,p,j,l,1));
                    end
                    intensity(num_ghost+n,num_ghost+p,j,l) = intensity(num_ghost+n,num_ghost+p,j,l) ...
                        -net_flux(num_ghost+n,num_ghost+p,j,l) -angular_flux(n,p,j,l); 
                end
            end
            %angular_flux = angular_flux + dt*v2/dphi*(circshift(i_flux(:,:,2),[0,-1]) - i_flux(:,:,2));

        end
     end %end of spatial loops
       
   
     
    %------------------------ NON-TIME SERIES OUTPUT ------------------ %
    if ~mod(i,output_interval)
        time
       % static_output(); 
        figure(2);            
        time_title = sprintf('t = %.3f (s)',time);
        %dont know where to add this above
        for j=1:nxa(1)
            %Ray intensity plots
            l=2; %select phi bin
            hi = subplot(2,3,j); 
            %since we pass matrices for the coordinates, do not transpose
            %intensity matrix
            h = pcolor(xx*ones(1,nx_r),ones(ny_r,1)*yy'+A*sin(xx*ones(1,nx_r).*K),intensity(is:ie,js:je,j,l));
            %turn off grid lines
            set(h, 'EdgeColor', 'none');
            
            %this string formatter doesnt work for some reason
            %subtitle = sprintf('$$\hat{k}^i_{Cartesian}$$ =(%0.3f, %0.3f,%0.3f)',mu(j,l,1),mu(j,l,2),mu(j,l,3));
            subtitle = sprintf('mu =(%0.3f, %0.3f,%0.3f)',mu(j,l,1),mu(j,l,2),mu(j,l,3));
            title(subtitle);
            %subtitle = ['$$\hat{k}^i_{Cartesian} = $$ (',num2str(mu(j,l,1),'%.3f'),',',num2str(mu(j,l,2),'%.3f'),',',...
            %    num2str(mu(j,l,3),'%.3f'),')'];
            %title(subtitle,'Interpreter','latex');      
            xlabel('x');
            ylabel('y + A sin(kx)');
            colorbar
        end
        %Mean intensity plot
         %figure(3);
         %pcolor(xx,yy,rad_energy(num_ghost+1:nx_r+2,num_ghost+1:ny_r+2)')
         pause(0.1)
    end
end
