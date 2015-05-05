%snake_angles_analytic.m
%Mini app that takes the formal solution to 
%------------------------ PARAMETERS ------------------ %
clear all;
%artificial division:
lx = 7*pi./(2*0.1); %k=0.1 here
ly = lx; 
nx = 200;
ny = nx;

c = 1.0;
dx = lx/nx;
dy = ly/ny;

N=12;
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
[ncells,nxa,mu,mu_b,pw] = uniform_angles2D(N,pi/2);
v = zeros(nx,ny,2); 
phi_bin = 1; %used to select phi level for plotting, injection of IC

%There is a notational confusion whereby nxa(2) includes the two poles.
%Therefore, there are actually nxa(2)-1 cells in the phi direction
%whereas in theta, we dont duplicate the last boundary ray (it is implied that 
% it is cyclic). Maybe want to be consistent in future
if (nxa(2) -1) >= 2
    num_phi_cells = nxa(2)-1; 
else
    num_phi_cells = 1;
end
%Radiation Spatial Discretization
%Radiation samples are volume averaged quantites
xx=linspace(0,lx,nx_r)';
yy=linspace(0,ly,ny_r)';

%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity_analytic = zeros(nx,ny,nxa(1),num_phi_cells); 

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

%------------------------ FORMAL SOLUTION ------------------------------ %

%Boundary conditions
%-------------------
for j=1:num_ghost
    %implement fixed dirichlet absorbing boundary conditions
    intensity_analytic(:,js-j,:) = 0.0;
    intensity_analytic(:,je+j,:) = 0.0;
    intensity_analytic(is-j,:,:) = 0.0;
    intensity_analytic(ie+j,:,:) = 0.0;
end
    
%Inject a ray in the -x +y direction from x_max, y_min
injection_ix = ie;
injection_jy = js +ny/2; 
injection_theta = 6; 
injection_phi = phi_bin;
for j=1:num_ghost
    intensity_analytic(injection_ix,injection_jy,injection_theta,injection_phi) = 1.0;
end
%technically long characteristics....
i = injection_ix-num_ghost;
j = injection_jy-num_ghost;
%Fill all cells in the path of the characteristic
x_pos = xx(i);
y_pos = yy(j);
x_dir = sign(mu(injection_theta,injection_phi,1)); %moving back or forward in x
y_dir = sign(mu(injection_theta,injection_phi,2));
slope_xy = mu(injection_theta,injection_phi,2)/mu(injection_theta,injection_phi,1); 
neighbor_x = i + x_dir;
neighbor_y = j + y_dir;

while (neighbor_x >= 1 && neighbor_y >= 1 && neighbor_x <= nx_r && neighbor_y <= ny_r) 
    %at one of four boundaries of a cell
    if x_dir > 0
        Dx = (xx(neighbor_x) - dx/2 - x_pos); 
    else
        Dx = (x_pos - xx(neighbor_x) - dx/2); 
    end
    
    if y_dir > 0
        Dy = (yy(neighbor_y) - dy/2 - y_pos);         
    else
        Dy = (y_pos - dy/2 - yy(neighbor_y));                
    end
    %see if characteristic hits x or y boundary first
    %x first
    if Dx/abs(mu(injection_theta,injection_phi,1)) <  ...
            Dy/abs(mu(injection_theta,injection_phi,2))
        i = neighbor_x; 
        neighbor_x = neighbor_x + x_dir;
        x_pos = xx(i);
        y_pos = y_pos + y_dir*abs(slope_xy)*Dx;
        intensity_analytic(i,j,injection_theta,injection_phi) = 1.0;
    else
        j = neighbor_y; 
        neighbor_y = neighbor_y + y_dir;
        y_pos = yy(j);
        x_pos = x_pos + x_dir*abs(Dy/slope_xy);
        intensity_analytic(i,j,injection_theta,injection_phi) = 1.0;
    end
    %CANNOT HANDLE CORNER ADVECTION
    
end

%------------------------ NON-TIME SERIES OUTPUT ------------------ %

h = figure(6);
clf;
set(h,'name','Analytic solution','numbertitle','off');
for j=1:nxa(1)
    %Ray intensity plots
    l=phi_bin; %select phi bin
    hi = subplot(3,4,j); 
    h = pcolor(xx,yy,intensity_analytic(num_ghost+1:nx_r+2,num_ghost+1:ny_r+2,j,l)');
    %turn off grid lines
    set(h, 'EdgeColor', 'none');
    %subtitle = ['$$\hat{k}^i_{Cartesian} = $$ (',num2str(mu(j,l,1),'%.3f'),',',num2str(mu(j,l,2),'%.3f'),',',...
    %    num2str(mu(j,l,3),'%.3f'),')'];
    %title(subtitle,'Interpreter','latex');      
    xlabel('x');
    ylabel('y');
    colorbar
end
%Mean intensity plot
 %figure(2);
 %pcolor(xx,yy,rad_energy(num_ghost+1:nx_r+2,num_ghost+1:ny_r+2)')

%TRANSFORM INTENSITY SOLUTION to snake angular bins (spatial cells should be the
%same)
    
%Renormalize the angular parameterization for snake as a function of
%spatial position
mu_s = zeros(nx_r,ny_r,nxa(1),num_phi_cells,3);
mu_b_s = zeros(nx_r,ny_r,nxa(1),nxa(2),3);
for i=1:nx_r
    beta_temp = beta(i);
    for j=1:ny_r
        for k=1:nxa(1)
            for l=1:num_phi_cells
                normalization_s = 1./sqrt((1+beta_temp.^2)*mu(k,l,1).^2 - 2*beta_temp*mu(k,l,1)*mu(k,l,2) + mu(k,l,2).^2 + mu(k,l,3).^2);
                normalization_b = 1./sqrt((1+beta_temp.^2)*mu_b(k,l,1).^2 - 2*beta_temp*mu_b(k,l,1)*mu_b(k,l,2) + mu_b(k,l,2).^2 + mu_b(k,l,3).^2);
                mu_s(i,j,k,l,:) = mu(k,l,:).*normalization_s; 
                mu_b_s(i,j,k,l,:) = mu_b(k,l,:).*normalization_b; 
                %check that each ray has unit spacelike norm in snake coords
                assert(abs(snake_norm(squeeze(mu_s(i,j,k,l,:)),sqrt(1+beta_temp^2),beta_temp) - 1.0) < normalization_tol);
                assert(abs(snake_norm(squeeze(mu_b_s(i,j,k,l,:)),sqrt(1+beta_temp^2),beta_temp)-1.0) < normalization_tol);
            end
            %Handle extra boundary line
            l=nxa(2);
            normalization_b = 1./sqrt((1+beta_temp.^2)*mu_b(k,l,1).^2 - 2*beta_temp*mu_b(k,l,1)*mu_b(k,l,2) + mu_b(k,l,2).^2 + mu_b(k,l,3).^2);
            mu_b_s(i,j,k,l,:) = mu_b(k,l,:).*normalization_b; 
            assert(abs(snake_norm(squeeze(mu_b_s(i,j,k,l,:)),sqrt(1+beta_temp^2),beta_temp) - 1.0) < normalization_tol )
        end
    end
end

intensity_analytic_s = zeros(nx_r,ny_r,nxa(1),num_phi_cells); 
for n=1:nx_r %should use overlap of solid angle bins
    for p=1:ny_r
        for j=1:nxa(1) % even though we only have one nonzero theta bin in current 
            %cartesian problem
            for l=num_phi_cells
                %Transform from Cartesian basis to Snake Basis
                mu_prime = zeros(3,1); 
                mu_prime(1) = mu(j,l,1);
                mu_prime(2) = (mu(j,l,2) + A*K*cos(K*xx(n))*mu(j,l,1));
                mu_prime(3) = mu(j,l,3);
                %check that transformed vector has norm 1 in snake coords
                assert(abs(snake_norm(mu_prime,sqrt(1+beta(n)^2),beta(n)) - 1.0) < normalization_tol);

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
                        intensity_analytic_s(n,p,k,l) = intensity_analytic_s(n,p,k,l) +...
                            intensity_analytic(n+num_ghost,p+num_ghost,j,l);
                      end
                end                         
            end
        end
    end
end

h = figure(7);
clf;
set(h,'name','Snake transformation of analytic solution','numbertitle','off');
for j=1:nxa(1)
    %Ray intensity plots
    l=phi_bin; %select phi bin
    hi = subplot(3,4,j); 
    h = pcolor(xx*ones(1,nx_r),ones(ny_r,1)*yy'+A*sin(xx*ones(1,nx_r).*K),intensity_analytic_s(:,:,j,l));
    %turn off grid lines
    set(h, 'EdgeColor', 'none');
    %subtitle = ['$$\hat{k}^i_{Cartesian} = $$ (',num2str(mu(j,l,1),'%.3f'),',',num2str(mu(j,l,2),'%.3f'),',',...
    %    num2str(mu(j,l,3),'%.3f'),')'];
    %title(subtitle,'Interpreter','latex');      
    xlabel('x');
    ylabel('y + A sin(kx)');
    colorbar
 end
 