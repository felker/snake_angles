function [num_rays,direction_cosines,mu_ordered,point_weights,level_weights] = angular_quad2D_snake(N)
%Input:
% N: Order of the quadrature = # of polar angles [0,pi] 
% only defined for N <= 12
%Output:
% num_rays
% direction_cosines 
% point_weights
% level_weights
%References: B.G. Carlson 1963, Bruls 1999
format long; 
%For 2D, we assume the physical domain is homogeneous in the z-axis.
%Therefore, we only do the positive hemisphere about the midplane, where
%the z-axis is the polar axis. 
%In 2D,
num_rays = N*(N+2)/2; %Up to 84 in 2d
num_rays_per_octant = num_rays/4;
direction_cosines = zeros(num_rays,2);
point_weights = zeros(num_rays,1);
level_weights = zeros(N/2,1); 
%This is the subset of [-1,1] that the direction_cosines are taken from in each
%dimension. Symmetry in each octant => only need one array
%Antisymmetric about zero
ordinates = zeros(1,N/2); 

%Choice of first cosine is arbitary. In analogy to Gaussian quadrature:
ordinates(1) = sqrt(1/(3*(N-1))); 
%In 3D space,
delta_cos =(1-3*ordinates(1)^2)* 2/(N-2); 

for i=2:N/2
    ordinates(i) = sqrt(ordinates(1)^2 + (i-1)*delta_cos);
end

%Derive weights for integration. See Bruls99 appendix
W = sqrt( 4/(3*(N-1)) + ([1:N/2] - 1)*2/(N-1));
%Compute level weights
level_weights(1) = W(1);
for i=2:N/2-1
    level_weights(i) = W(i) - W(i-1);
end

%following ATHENA, we correct the last level weight to renormalize...not
%sure why this happens
level_weights(N/2) = 1.0 - sum(level_weights); 
%permuation matrix
pl = zeros(N/2,3);
pmat = zeros(N/2,N/2); 
plab = zeros(num_rays,1);

%To select valid rays for the quadrature, we follow the index selection
%rule:
%sigma = i + j + k, indices of the direction cosines must equal ntheata/2 +2
%for normalization in 3D space
ray_count = 0;
np = 0;
ip = 0;
for i=1:N/2
   for j=1:(N/2+1-i)
       k = N/2+2-i-j;
       %how do we order the rays in each octant?
       ray_count = ray_count +1;
       direction_cosines(ray_count,1) = ordinates(i);
       direction_cosines(ray_count,2) = ordinates(j);
       if N <= 12 %system of equations is lin. ind. iff N<=12
       ip = permutation(i,j,k,pl,np);
       if ip == -1 %ray indices havent been loaded yet
           np = np+1;
           pl(np,1) = i;
           pl(np,2) = j;
           pl(np,3) = k;
           pmat(i,np) = pmat(i,np) + 1.0;
           plab(ray_count) = np;
       else %these indices have been loaded 
          pmat(i,ip) = pmat(i,ip) + 1.0;
          plab(ray_count) = ip; 
       end
       end
   end
end
assert(ray_count == num_rays_per_octant);

%Symmetry: reflect across second axis 
%(must be a better way of doing this)
direction_cosines(ray_count+1:2*ray_count,1) = -direction_cosines(1:ray_count,1);
direction_cosines(ray_count+1:2*ray_count,2) = direction_cosines(1:ray_count,2);
ray_count = ray_count*2;
%Symmetry: reflect across first axis
direction_cosines(ray_count+1:2*ray_count,1) = direction_cosines(1:ray_count,1);
direction_cosines(ray_count+1:2*ray_count,2) = -direction_cosines(1:ray_count,2);
ray_count = ray_count*2;
%Material symmetry: do not reflect across final axis

%test ray normalization
%some of these norms have a deviation from 1.0 in the 15th decimal place
%assert(all((abs(sum((direction_cosines.*direction_cosines),2) - ones(num_rays,1)) < 1e-12)));

%test output in unit sphere
%quiver3(zeros(ray_count,1),zeros(ray_count,1),zeros(ray_count,1), ...
 %   direction_cosines(1:ray_count,1),direction_cosines(1:ray_count,2), ...
  %  direction_cosines(1:ray_count,3)); 
%test projection onto z=0 plane
%quiver(zeros(num_rays,1),zeros(num_rays,1),direction_cosines(:,1),direction_cosines(:,2)); 

%Use Bruls method to calculate families of weights. Solve system of eqs
wpf = pmat(1:N/2-1,1:N/2-1)\level_weights(1:N/2-1);
for i=1:num_rays_per_octant %each quadrant, normalize %Is this ordered correctly?
    point_weights(i) = wpf(plab(i))*0.25;
    point_weights(i+num_rays_per_octant) = wpf(plab(i))*0.25; 
    point_weights(i+2*num_rays_per_octant) = wpf(plab(i))*0.25; 
    point_weights(i+3*num_rays_per_octant) = wpf(plab(i))*0.25; 
end

%ADDED FOR SNAKE SOLVER
%logical array to find neighboring solid angles
mu_ordered = NaN(N,N,2); 
    for i=1:N/2 %fill one half of the jagged array
        for j=1:i %one quadrant is a lower triangular matrix 
            k = N/2 + j; 
            l = N/2 + 1 -j;
            mu_ordered(i,k,:) = [ordinates(j), ordinates(N/2 + 1 - i)];
            mu_ordered(i,l,:) = [-ordinates(j), ordinates(N/2 +1 - i)];
        end
    end
    A = flipud(mu_ordered(:,:,1));
    B = -flipud(mu_ordered(:,:,2));
    mu_ordered(N/2+1:N,:,1) = A(N/2+1:N,:);
    mu_ordered(N/2+1:N,:,2) = B(N/2+1:N,:);
%test output in unit sphere at two example points in snake system
close all
%one at x',y' = pi/2k, A*sin(kx) = pi/2k, A
%d_x basis vector = (1,-Akcos(kx)) = (1,0)
%----> SAME AS CARTESIAN GLOBAL TANGENT SPACE

mu_z = sqrt(1- direction_cosines(1:ray_count,1).^2 - direction_cosines(1:ray_count,2).^2);
quiver3(zeros(ray_count,1),zeros(ray_count,1),zeros(ray_count,1), ...
   direction_cosines(1:ray_count,1),direction_cosines(1:ray_count,2), ...
   mu_z,0);
%draw level circles
hold on;
xlabel('x_{Cartesian}');
ylabel('y_{Cartesian}');
zlabel('z_{Cartesian}');
quiver3(0,0,0,1,0,0,0,'-r');

%Now, we plot the level circles ala Carlsson. This is treating the rays as
%nodal values
% t=0:0.01:2*pi;
% st = sin(t);
% ct = cos(t);
% for i=1:ray_count %in x-y plane
%     radius = sqrt(direction_cosines(i,1).^2 + direction_cosines(i,2).^2);
%     plot3(radius*st,radius*ct,mu_z(i)*ones(size(t)),'-k'); 
% end
% t=0:0.01:pi;
% st = sin(t);
% ct = cos(t);
% for i=1:ray_count %in x-z plane
%     radius = sqrt(direction_cosines(i,1).^2 + mu_z(i).^2);
%     plot3(radius*ct,direction_cosines(i,2)*ones(size(t)),radius*st,'-k'); 
% end
% for i=1:ray_count %in y-z plane
%     radius = sqrt(direction_cosines(i,2).^2 + mu_z(i).^2);
%     plot3(direction_cosines(i,1)*ones(size(t)),radius*ct,radius*st,'-k'); 
% end

% If we interpret each ray as a volume averaged solid angle quantity, how
% can we visualize the cells?
%Start with a meshing specific to Carlsson discretization. %Put circles at
%angles halfway between the current angles
t=0:0.01:pi;
st = sin(t);
ct = cos(t);
%ONLY NEED TO DO ONE HALF OF MATRIX
% for i=1:N/2 %loop down rows, difference from left to right Each row has same mu_y
%     for j=1:i %go from center to right
%         if (i ~= 1 && j==i) %after the first row, have to loop around
%             k = N/2 + j; 
%             mu_x_mid = cos((acos(mu_ordered(i,k,1)) + acos(mu_ordered(i,N+1-i,1)))/2);
%         else
%             k = N/2 + j; 
%             l= k-1;
%             %plot in y-z plane
%             mu_x_mid = cos((acos(mu_ordered(i,k,1)) + acos(mu_ordered(i,l,1)))/2);
%         end
%         radius = sqrt(1-mu_x_mid^2); 
%         plot3(mu_x_mid*ones(size(t)),radius*ct,radius*st,'-k')
%     end
%     for m=1:i %go from center to left  
%             k = N/2 + 1 -m; 
%             l= k+1;
%             %plot in y-z plane
%             mu_x_mid = cos((acos(mu_ordered(i,k,1)) + acos(mu_ordered(i,l,1)))/2)
%                     radius = sqrt(1-mu_x_mid^2); 
%         plot3(mu_x_mid*ones(size(t)),radius*ct,radius*st,'-k')
%     end
% end
% for i=1:N/2 %loop across columns, difference from top to bottom. Each col has same mu_x
%     for j=1:i %go from center to top
%         if (i ~= 1 && j==i) %after the first column, have to loop around
%             k = N/2 + j; 
%             mu_y_mid = cos((acos(mu_ordered(k,i,2)) + acos(mu_ordered(N+1-i,k,2)))/2);
%         else
%             k = N/2 + j; 
%             l= k-1;
%             %plot in y-z plane
%             mu_y_mid = cos((acos(mu_ordered(k,i,2)) + acos(mu_ordered(l,i,2)))/2);
%         end
%         radius = sqrt(1-mu_y_mid^2); 
%         plot3(radius*ct,mu_y_mid*ones(size(t)),radius*st,'-k')
%     end
%     for m=1:i %go from 
%             k = N/2 + 1 -m; 
%             l= k+1;
%             %plot in x-z plane
%             mu_y_mid = cos((acos(mu_ordered(k,i,2)) + acos(mu_ordered(l,i,2)))/2)
%                     radius = sqrt(1-mu_y_mid^2); 
%         plot3(radius*ct,mu_y_mid*ones(size(t)),radius*st,'-k')
%     end
% end

%Other representative extreme at x',y' = 0,0
%d_x basis vector = (1,-Akcos(kx)) = (1,-Ak)
%ADD -Ak to y-component of sampling angles?
%need to renormalize! 
% A=10.0;
% k=0.1;
% hold off;
% figure(2);
% %this C normalization constant is a vector of constants
% normalization = 1.0./sqrt(1+A^2*k^2.*direction_cosines(1:ray_count,1).^2 - 2*A*k*...
%     direction_cosines(1:ray_count,1).*direction_cosines(1:ray_count,2));
% direction_cosines(1:ray_count,1) = normalization.*(direction_cosines(1:ray_count,1));
% direction_cosines(1:ray_count,2) = normalization.*(direction_cosines(1:ray_count,2) -A*k*direction_cosines(1:ray_count,1));
% 
% mu_z = normalization.*mu_z;
% %these vectors are  normalized in cartesian coordinates, but are NOT
% %normalized in snake metric
% quiver3(zeros(ray_count,1),zeros(ray_count,1),zeros(ray_count,1), ...
%    direction_cosines(1:ray_count,1),direction_cosines(1:ray_count,2), ...
%    mu_z,0);
% xlabel('x_{Cartesian}');
% ylabel('y_{Cartesian}');
% zlabel('z_{Cartesian}');
% hold on;
% quiver3(0,0,0,1,0,0,0,'-r');
% quiver3(0,0,0,1./sqrt(1+A^2*k^2),-A*k./sqrt(1+A^2*k^2),0,0,'-r');

%draw level circles
% t=0:0.01:2*pi;
% st = sin(t);
% ct = cos(t);
% for i=1:ray_count %in x-y plane
%     radius = sqrt(direction_cosines(i,1).^2 + direction_cosines(i,2).^2);
%     plot3(radius*st,radius*ct,mu_z(i)*ones(size(t)),'-k'); 
% end
% t=0:0.01:pi;
% st = sin(t);
% ct = cos(t);
% for i=1:ray_count %in x-z plane
%     radius = sqrt(direction_cosines(i,1).^2 + mu_z(i).^2);
%     plot3(radius*ct,direction_cosines(i,2)*ones(size(t)),radius*st,'-k'); 
% end
% for i=1:ray_count %in y-z plane
%     radius = sqrt(direction_cosines(i,2).^2 + mu_z(i).^2);
%     plot3(direction_cosines(i,1)*ones(size(t)),radius*ct,radius*st,'-k'); 
% end


%THIS STUFF IS TAKEN FROM SNAKE_ANGLES.M FOR CALUCLATING FLUXES FROM
%MU_ORDERED-- WHICH WE ARE NO LONGER USING
    %Hardcode the donor cell method for solid angle cells at each spatial
    %point
%     for n=1:nx_r 
%         for p=1:ny_r
%             %Longitudanal flux (dimension 1 in array)
%             for j=1:ntheta/2 %positive mu_y
%                 for k=1:j %one quadrant is a lower triangular matrix
%                     l = ntheta/2 + k; %positive mu_x
%                     m = ntheta/2 + 1 -k; %negative mu_x
%                     if k==j %we have reached the ends of row of mu_ordered
%                     end
%                     %Need to find the lattitudanal neighbors at the ends
%                     %of the array?
%                 end
%             end
%                     
%             %Lattitudanal flux (dimension 2 in array)
%             %handle coordinate singularity at poles 
%         end
%     end
%            for j=1:ntheta/2 %positive mu_y
%                 for k=1:j %one quadrant is a lower triangular matrix 
%                     l = ntheta/2 + k; %positive mu_x
%                     m = ntheta/2 + 1 -k; %negative mu_x
%                     C2(n,p,j,l) = A*k^2*sin(k*xx(n))*(1+2*beta(n).^2)*mu_ordered(j,l,1).^2; %beta(x)
%                     C2(n,p,j,m) = A*k^2*sin(k*xx(n))*(1+2*beta(n).^2)*mu_ordered(j,m,1).^2;
%                 end
%             end
%             for j=ntheta:-1:ntheta/2+1 %negative mu_y
%                 for k=1:(ntheta-j+1) %one quadrant is a lower triangular matrix 
%                     l = ntheta/2 + k; %positive mu_x
%                     m = ntheta/2 + 1 -k; %negative mu_x
%                     %compute C^2
%                     C2(n,p,j,l) = A*k^2*sin(k*xx(n))*(1+2*beta(n).^2)*mu_ordered(j,l,1).^2; %beta(x)
%                     C2(n,p,j,m) = A*k^2*sin(k*xx(n))*(1+2*beta(n).^2)*mu_ordered(j,m,1).^2;
%                 end
%             end
end