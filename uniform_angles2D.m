function [ncells,nxa,mu,mu_b,pw] = uniform_angles2D(N)
%UNIFORM_ANGLES2D Generate uniform 2D discretization of mu
%   Specifically created for the Snake Coords solver, this simply returns
%   a mesh of the 2-sphere mu_x,mu_y,mu_z \in [-1,1] x [-1,1] x [0,1]
%   with the constraint that mu_x^2 + mu_y^2 + mu_z^2 = 1 
%
%   Input: N-- Number of mu level cells in one \hat{k}^i dimension (N/2 for
%   mu_z in 2D) (note, there are n+1, n/2+1 boundaries, including -1,1 and
%   0,1) N MUST BE EVEN, shoul be greater than 6

%   We can specify each Nkx,nky,nkz separately. For now, assume N= nkx=nky
% then, each level will have 2*(N-1). And there are N/2 - 1 levels
% therefore, ncells=2*(N-1)*(N/2 -1)

%   Since a uniform discretization will cluster rays near poles, best
%   to use with large N for good accuracy (but inefficient)
%
%   Output: ncells: total number of FV solid angle clells 
%   nxa: angular parameter coordinate dimensions. currently, it is 2-dim
%   with nxa(1) = N nxa(2) = N/2 +1.
%   

% dk = 2.0./(N); %increments in the z-coordinate. Must have a cell terminate
% %on the pole (mu_z=1!)
% mu_z = [0:dk:1];
% 
% %produces N/2-1 z cells + 1 polar cap NOT ANYMORE dk3 = 2.0./(N+1)
% 
% %Now, given the z-levels, discretize that strip into a ring of solid angle
% %cells
% 
% %we make this easy by having a cell boundary at mu_z=0. This allows N to
% %specify the number of cells in a unit circle. Then all higher mu_z levels
% %are fixed by this initial discretization
% close all;
% for i=1:(N/2) %foreach mu_z level
%     remainder = 1 - mu_z(i).^2;
%     mu_x = [-remainder:dk:0];
%     %second and third quadrants
%     mu_y =  [sqrt(remainder - mu_x.^2), -sqrt(remainder - mu_x.^2) ];
%     mu_x = [mu_x, [0:dk:remainder]]; 
%     mu_y = [mu_y, sqrt(remainder - mu_x(N/2+1:N+2).^2), -sqrt(remainder - mu_x(N/2+1:N+2).^2) ];
%     % mu_y =  [sqrt(remainder - mu_x(1:N/2+1).^2),-sqrt(remainder - fliplr(mu_x(1:N/2+1).^2))] 
%     %N+2 levels due to repetition of 0
%     t=0:0.01:2*pi;
%     st = sin(t);
%     ct = cos(t);
%     hold on;
%     plot(ct,st,'-k');
%     q= quiver(zeros(1,N+2),zeros(1,N+2),mu_x,mu_y,0,'.');
% end
% hold off; 

%NEVERMIND--- WE WANT TO DISCRETIZE THE PARAMETERS: THETA, PHI
%In 2D problems,theta \in [0, pi], \phi
%Lets just do the entire 2-sphere

% theta = linspace(0,2.0*pi,N);
% phi = linspace(0,pi,N/2);
%dont mesh all the way to 2pi

%now, keep discretization equal for ease of sphere(n)
%nxa = [N, N/2 +1];
nxa = [N, N+1];
dtheta = 2*pi./(nxa(1)); 
dphi = pi./(nxa(2)-1);
%dphi = pi./(N/2); 
theta= [0:dtheta:2.0*pi-dtheta];
phi= [0:dphi:pi];

%rows are theta bin, columns are phi bin
kx = cos(theta)'*sin(phi);
ky = sin(theta)'*sin(phi);
kz = ones(N,1)*cos(phi);

%plot circles at each mu_z in x-y plane
t=0:0.01:2*pi;
st = sin(t);
ct = cos(t);
figure(1);
hold on;
for i=1:nxa(2) %at fixed phi. in x-y plane, doesnt include poles
     radius = sqrt(1-cos(phi(i)).^2);
     plot3(radius*st,radius*ct,cos(phi(i))*ones(size(t)),'-k'); 
end
%these semicircles are at fixed theta. 
t=0:0.01:pi;
st = sin(t);
ct = cos(t);
for i=1:nxa(1)
     radius = 1.0;
     plot3(cos(theta(i))*st,sin(theta(i))*st,radius*ct,'-k'); 
end

for i=1:nxa(1)
    for j=1:nxa(2)
        quiver3(0,0,0,kx(i,j),ky(i,j),kz(i,j),0,'.b')
    end
end

%periodicity of the domain implies duplication of 1 theta bin
%first and last column/row are identical for kx, ky
%how do we double the number of cells in the kx,ky,kz structures?
%ncells=(N-1)*(N/2 -1);
%pw = 1./((N-1)*(N/2-1))*ones(2*(N-1),N/2-1);

ncells= nxa(1)*(nxa(2)-1);
mu_b = zeros(nxa(1),nxa(2),3);
%these are the \hat{k} boundaries of each solid angle cell

%THESE ARE IN SNAKE BASIS, NOT CARTESIAN!!
mu_b(:,:,1) = kx;
mu_b(:,:,2) = ky;
mu_b(:,:,3) = kz;
%we need the sampling angle for each FV cell, must be unit spacelike
%vectors
%lets use average of angles to ensure normalizaiton
mu = zeros(nxa(1),nxa(2)-1,3); %identify with left, top boundaries
for i=1:nxa(1)-1
    for j=1:nxa(2)-1
        mu(i,j,1) = cos((theta(i) + theta(i+1))./2)*sin((phi(j) + phi(j+1))./2);
        mu(i,j,2) = sin((theta(i) + theta(i+1))./2)*sin((phi(j) + phi(j+1))./2);
        mu(i,j,3) = cos((phi(j) + phi(j+1))./2);
    end
end
%wrap around theta dimension. the first angle should be theta=2pi, not
%theta=0. This matters for averaging, otherwise you wont get hte last theta
%ray
for j=1:nxa(2)-1
        mu(N,j,1) = cos((theta(N) + 2*pi)./2)*sin((phi(j) + phi(j+1))./2);
        mu(N,j,2) = sin((theta(N) + 2*pi)./2)*sin((phi(j) + phi(j+1))./2);
        mu(N,j,3) = cos((phi(j) + phi(j+1))./2);
end
%plot representative vectors
for i=1:nxa(1)
    for j=1:nxa(2)-1
        quiver3(0,0,0,mu(i,j,1),mu(i,j,2),mu(i,j,3),0,'r')
    end
end
xlabel('x');
ylabel('y');
zlabel('z');
hold off;
%quadrature weights. Should this be proportional to the "size of the cell"?
pw = 1./(ncells)*ones(nxa(1),nxa(2)-1);
%the trick is ensuring unit norm of each cell boundary... if we were
%discretizing the unit cube, we would independently discretize each
%direction, since the coordinates are wholly independent. THIS IS WHY
%DISCRETIING THE 2 DOF PARAMETERS IS SIMPLER

%for the 2-sphere, we fix the discretizaiton of one dimension, and then the
%other two are determined by the remainder on the circle.

end

