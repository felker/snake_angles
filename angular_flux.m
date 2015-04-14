%angular_flux.m
%test script to solve the advection equation in solid angle
%du/dt + v1 du/dxA(1) + v2 du/dxA(2) = 0

%this will help me put the full solver together

N=8;
nt = 100;
dt = 0.1;

[ncells,nxa,mu,mu_b,pw] = uniform_angles2D(N);

u = zeros(nxa(1),nxa(2)); 

%initial condition
u(1,4) = 1.0; % initially = away from poles
v1=0.1;
v2=0.0;

%what is the CFL condition for the angular bins?
dtheta = 2*pi./(nxa(1)); 
dphi = pi./(nxa(2));

cfl_angular = dt*[v1 v2]* [1./dtheta 1./dphi]';
dt = dt*1.0/cfl_angular;

hold off; 
for i=1:nt 
   time =i*dt %#ok<NOPTS>
   %plot solid angle surface volume plot
   figure(2);
   [x,y,z] = sphere(N);
   surf(x,y,z,u');
   colorbar;
   caxis([-1.0 1.0])

   net_flux = zeros(nxa(1),nxa(2));
   i_flux = zeros(nxa(1),nxa(2),2);
        %theta flux
        if v1 > 0 %counterclockwise speed
            i_flux(:,:,1) = circshift(u(:,:),[1,0]); %flux entering left boundary
        elseif v1 < 0 %clockwise 
            i_flux(:,:,1) = u(:,:);  %flux leaving left boundary          
        end
        %phi flux
        for j=2:nxa(2)-1
            if v2 > 0 
                i_flux(:,j,2) = circshift(u(:,j),[0,1]); 
            elseif v2 < 0
                i_flux(:,j,2) = u(:,j);            
            end
        end
        %handle both poles
%         j=1; 
%         if v2 > 0
%         j=nxa(2);
        
   net_flux = dt*v1/dtheta*(circshift(i_flux(:,:,1),[1,0]) - i_flux(:,:,1));
   net_flux = net_flux + dt*v2/dphi*(circshift(i_flux(:,:,2),[0,1]) - i_flux(:,:,2));
   u = u - net_flux; 
   pause();
end