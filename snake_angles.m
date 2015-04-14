%snake_angles.m
%Mini-app that does vacuum RTE. we inject an initial propagation of rays
%and then transform these rays to snake coordinates to see how the angular
%mesh evolves. 

%this only tests advection step

%------------------------ PARAMETERS ------------------ %
clear all;
close all;
nx = 32;
ny = 32;
N = 4; %quadarture is only defined up to 12 in each direction, must be even
%order of the quadrature, N, refers to the number of mu-levels in the interval [-1, 1].

lx = 100.0;
ly = 100.0; %corresponding to max y_cartesian, not the actual max y' 
c = 1.0;
dx = lx/nx;
dy = ly/ny;
dt = 0.0025;
nt = 32;
%Upwind monotonic interpolation scheme
method = 'van Leer'; 
time_centering = 0; %explicit = 0, implicit = 1, CN=1/2

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
v = zeros(nx,ny,2); 

%Radiation Spatial Discretization
%Radiation samples are volume averaged quantites
xx=linspace(0,lx,nx_r)';
yy=linspace(0,ly,ny_r)';

%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity = zeros(nx,ny,nxa(1),nxa(2)); 

%-------------------- SNAKE COORDINATE DETAILS  ----------------- %
A = 10.0; %amplitude of snake oscillations in y'
k = 0.1; %wavenumber of oscillations 
beta = A*k*cos(k*xx);
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
%THIS NEEDS TO BE CALCULATED AT EACH POINT!
cfl_mu = C*dt*abs(mu)*[1/dx 1/dy]';
%set dt to have max cfl of 0.4
dt = dt*1.0/max(cfl_mu);
%Recalculate radiation CFL
cfl_mu = C*dt*abs(mu)*[1/dx 1/dy]';
assert(min(abs(cfl_mu) <= ones(na,1))); 

dt = 0.0025;

%------------------------ VELOCITY TERMS ------------------------------ %
[nv, vvnn, vCsquare, vsquare, absV] = update_velocity_terms(v,mu,C);

%------------------------ OUTPUT VARIABLES------------------------------ %
output_interval = 2; 
num_output = 8; %number of data to output
num_pts = nt/output_interval; 
time_out = dt*linspace(0,nt+output_interval,num_pts+1); %extra pt for final step
y_out = zeros(num_pts+1,num_output);

%Explicit-Implicit operator splitting scheme
%-------------------------------------------
for i=0:nt
    time = dt*i %#ok<NOPTS> 
    %Substep 0: Time series output, boundary condition, moments update
    %Update moments
    [J,H,K,rad_energy,rad_flux,rad_pressure] = update_moments(intensity,mu,pw,c);
    photon_momentum = P*rad_flux/C;
    if ~mod(i,output_interval)
      %  [y_out] = time_series_output(y_out,time_out,rad_energy,rad_flux,rad_pressure,v,vCsquare,nx,ny,C,GasMomentum,GasKE,GasTE,photon_momentum);
    end

    %Inject a single ray in the  mu=(1/3,1/3) from center
    intensity(nx/2,ny/2,1) = 1.0; 
    
    %Box periodic boundary conditions
    for j=1:num_ghost
        intensity(:,js-j,1:na/2) = intensity(:,je+1-j,1:na/2); 
        intensity(:,je+j,na/2:na) = intensity(:,js-1+j,na/2:na); 
        %mu_x arent easily divided into positive and negative
        intensity(is-j,:,1:3) = intensity(ie+1-j,:,1:3); 
        intensity(is-j,:,7:9) = intensity(ie+1-j,:,7:9); 
        intensity(ie+j,:,4:6) = intensity(is-1+j,:,4:6); 
        intensity(ie+j,:,10:12) = intensity(is-1+j,:,10:12); 
    end
    
    %Substep #1: Explicitly advance transport term
    net_flux = zeros(nx,ny,na);
    for j=1:na %do all nx, ny at once
        %Split the transport term into diffusion and advection terms
        %Photon diffusion at nearly the speed of light
        tau = (10*(dx)*(rho_a+rho_s)).^2; %Eq 14 in Jiang14
        %Correct for division by zero 
        for k=1:nx
            for l=1:ny
                if tau(k,l) > 0 
                    alpha(k,l) = sqrt((1-exp(-tau(k,l)))./(tau(k,l)));
                else
                    alpha(k,l) = 1.0;
                end
            end
        end
       
        %Advection with the fluid only if there are material terms
        advection_term = (3*nv(:,:,j).*J).*(abs(rho_a+rho_s)>0);
        i_flux = upwind_interpolate2D_decoupledV2(mu(j,1)*(intensity(:,:,j) - advection_term/C),method,dt,dx,alpha*C*sign(mu(j,1)),is,ie+1,js,je+1,1);
        net_flux(is:ie,js:je,j) = dt*C/dx*(i_flux(is+1:ie+1,js:je) - i_flux(is:ie,js:je));
        %should also turn off advection velocity if the projection along
        %axis is zero in rare cases
        %Also, should take into account velocity gradient between cells
        advection_term = (mu(j,1)^2*3*J).*(abs(rho_a+rho_s)>0);
        velocity_term = v(:,:,1) + mu(j,2)*v(:,:,2)/mu(j,1);
        i_flux = upwind_interpolate2D_decoupledV2(advection_term, method, dt, dx, velocity_term,is,ie+1,js,je+1,1);
        net_flux(is:ie,js:je,j) = net_flux(is:ie,js:je,j) + ...
            dt/dx*0.5*((velocity_term(is+1:ie+1,js:je)+velocity_term(is:ie,js:je))*i_flux(is+1:ie+1,js:je) -...
            (velocity_term(is-1:ie-1,js:je)+velocity_term(is:ie,js:je))*i_flux(is:ie,js:je));
    end %end of ray trace loop, x-direction
    
    for j=1:na %do all nx, ny at once
        %Split the transport term into diffusion and advection terms
        %Photon diffusion at nearly the speed of light
        tau = (10*(dy)*(rho_a+rho_s)).^2; %Eq 14 in Jiang14
        %Correct for division by zero 
        for k=1:nx
            for l=1:ny
                if tau(k,l) > 0 
                    alpha(k,l) = sqrt((1-exp(-tau(k,l)))./(tau(k,l)));
                else
                    alpha(k,l) = 1.0;
                end
            end
        end
         %Advection with the fluid only if there are material terms
        advection_term = (3*nv(:,:,j).*J).*(abs(rho_a+rho_s)>0);
        i_flux = upwind_interpolate2D_decoupledV2(mu(j,2)*(intensity(:,:,j) - advection_term/C),method,dt,dy,alpha*C*sign(mu(j,2)),is,ie+1,js,je+1,2);
        net_flux(is:ie,js:je,j) = net_flux(is:ie,js:je,j)+ dt*C/dy*(i_flux(is:ie,js+1:je+1) - i_flux(is:ie,js:je));
        %should also turn off advection velocity if the projection along
        %axis is zero in rare cases
        %Also, should take into account velocity gradient between cells
        advection_term = (mu(j,2)^2*3*J).*(abs(rho_a+rho_s)>0);
        velocity_term = v(:,:,2) + mu(j,1)*v(:,:,1)/mu(j,2);
        i_flux = upwind_interpolate2D_decoupledV2(advection_term, method, dt, dy, velocity_term,is,ie+1,js,je+1,2);
        net_flux(is:ie,js:je,j) = net_flux(is:ie,js:je,j) + ...
            dt/dy*0.5*((velocity_term(is:ie,js+1:je+1)+velocity_term(is:ie,js:je))*i_flux(is:ie,js+1:je+1) -...
            (velocity_term(is:ie,js-1:je-1)+velocity_term(is:ie,js:je))*i_flux(is:ie,js:je));
    end    %end of ray trace loop, y-direction
    %Substep #1.1: Compute solid angular fluxes
    %the coefficients C^i can depend C^i(x,y,z,k^1,k^2,k^3), so the outer
    %loop should be over solid angle, and the inner loops will be over
    %spatial cells... actually the order shouldnt matter
        
    %Compute C^2
    %for this metric/coords, this can be precomputed. varies from -1C2 to
    %1C2
    C2 = NaN(nx_r,ny_r,ntheta,ntheta);
    for n=1:nx_r 
        for p=1:ny_r
        %split the loop into two sweeps of the mu_ordered matrix
            for j=1:ntheta/2 %positive mu_y
                for k=1:j %one quadrant is a lower triangular matrix 
                    l = ntheta/2 + k; %positive mu_x
                    m = ntheta/2 + 1 -k; %negative mu_x
                    C2(n,p,j,l) = A*k^2*sin(k*xx(n))*(1+2*beta(n).^2)*mu_ordered(j,l,1).^2; %beta(x)
                    C2(n,p,j,m) = A*k^2*sin(k*xx(n))*(1+2*beta(n).^2)*mu_ordered(j,m,1).^2;
                end
            end
            for j=ntheta:-1:ntheta/2+1 %negative mu_y
                for k=1:(ntheta-j+1) %one quadrant is a lower triangular matrix 
                    l = ntheta/2 + k; %positive mu_x
                    m = ntheta/2 + 1 -k; %negative mu_x
                    %compute C^2
                    C2(n,p,j,l) = A*k^2*sin(k*xx(n))*(1+2*beta(n).^2)*mu_ordered(j,l,1).^2; %beta(x)
                    C2(n,p,j,m) = A*k^2*sin(k*xx(n))*(1+2*beta(n).^2)*mu_ordered(j,m,1).^2;
                end
            end
        end
    end
    return;
    angular_flux = zeros(nx_r,ny_r,ntheta,ntheta,2);
    %Hardcode the donor cell method for solid angle cells at each spatial
    %point
    for n=1:nx_r 
        for p=1:ny_r
            %Longitudanal flux (dimension 1 in array)
            for j=1:ntheta/2 %positive mu_y
                for k=1:j %one quadrant is a lower triangular matrix
                    l = ntheta/2 + k; %positive mu_x
                    m = ntheta/2 + 1 -k; %negative mu_x
                    if k==j %we have reached the ends of row of mu_ordered
                    end
                    %Need to find the lattitudanal neighbors at the ends
                    %of the array?
                end
            end
                    
            %Lattitudanal flux (dimension 2 in array)
            %handle coordinate singularity at poles 
        end
    end

    intensity = intensity - net_flux; 
    [J,H,K,rad_energy,rad_flux,rad_pressure] = update_moments(intensity,mu,pw,c);
    %------------------------ NON-TIME SERIES OUTPUT ------------------ %
    if ~mod(i,output_interval)
       % static_output(); 
        for j=1:na/2
                   hi = subplot(2,3,j); 
                   h = pcolor(xx,yy,intensity(num_ghost+1:nx_r+2,num_ghost+1:ny_r+2,j)');
                   set(h, 'EdgeColor', 'none');
                   x_label = sprintf('mu =(%f, %f)',mu(j,1),mu(j,2));
                   time_title = sprintf('t = %f (s)',time);
                   xlabel(x_label);
                   title(time_title);
                   colorbar
                      end
                figure(2);
        pcolor(xx,yy,rad_energy(num_ghost+1:nx_r+2,num_ghost+1:ny_r+2)')
               pause(0.5)
    end
    
end
