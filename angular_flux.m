%angular_flux.m
%test script to solve the advection equation in solid angle
%du/dt + v1 du/dxA(1) + v2 du/dxA(2) = 0

%this will help me put the full solver together

N=8;
nt = 2*pi*1000; %needs to go 2pi
dt = 0.01;

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

hold off; 
for i=1:nt 
   time =i*dt %#ok<NOPTS>
   %plot solid angle surface volume plot
   if mod(i,10) == 0
       figure(2);
       [x,y,z] = sphere(N);
       surf(x,y,z,u');
       colorbar;
       caxis([-1.0 1.0])
       xlabel('x');
       ylabel('y');
       zlabel('z');
       M(i./10) = getframe;
   end
   
   net_flux = zeros(nxa(1),nxa(2));
   i_flux = zeros(nxa(1),nxa(2),2);
        %theta flux
        if v1 > 0 %counterclockwise speed
            i_flux(:,:,1) = circshift(u(:,:),[1,0]); %flux entering "left" boundary
        elseif v1 < 0 %clockwise 
            i_flux(:,:,1) = u(:,:);  %flux leaving left boundary          
        end
        %phi flux
        for j=2:nxa(2)-1
            if v2 > 0 
                i_flux(:,j,2) = u(:,j-1); 
            elseif v2 < 0
                i_flux(:,j,2) = u(:,j);            
            end
        end
        %handle both poles
        j=1;
        if v2 > 0
            for k=1:nxa(1)/2 %find opposite theta bin
                i_flux(k,j,2) = u(k,j);
            end
            for k=nxa(1)/2+1:nxa(1)
                i_flux(k,j,2) = u(mod(k+nxa(1)/2,nxa(1)),j);
            end
        elseif v2 < 0 
            for k=1:nxa(1)/2 %find opposite theta bin
                i_flux(k,j,2) = u(k+nxa(1)/2,j);
            end
            for k=nxa(1)/2+1:nxa(1)
                i_flux(k,j,2) = u(k,j);
            end
        end
        j=nxa(2);
        if v2 > 0
            for k=1:nxa(1)/2 
                i_flux(k,j,2) = u(k,j-1);
            end
            for k=nxa(1)/2+1:nxa(1)%find opposite theta bin
                i_flux(k,j,2) = u(mod(k+nxa(1)/2,nxa(1)),j);
            end
        elseif v2 < 0 
            for k=1:nxa(1)/2 
                i_flux(k,j,2) = u(k+nxa(1)/2,j);
            end
            for k=nxa(1)/2+1:nxa(1)
                i_flux(k,j,2) = u(k,j);
            end
        end        
   net_flux = dt*v1/dtheta*(circshift(i_flux(:,:,1),[-1,0]) - i_flux(:,:,1));
   net_flux = net_flux + dt*v2/dphi*(circshift(i_flux(:,:,2),[0,-1]) - i_flux(:,:,2));
   u = u - net_flux; 
   %pause();
end
pause();
figure
movie(M)