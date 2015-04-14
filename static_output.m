function [] = static_output()
%STATIC_OUTPUT Used to create plots that don't collect data at multiple
%timesteps. 

%     figure(2);
 %      pcolor(xx(2:nx-1,1),yy(2:ny-1,1),rad_energy(2:nx-1,2:ny-1)')
   %colorbar
%         hi = subplot(2,3,1); 
%         h = pcolor(xx,yy,intensity(:,:,1)');
%         set(h, 'EdgeColor', 'none');
%         x_label = sprintf('mu =(%f, %f)',mu(1,1),mu(1,2));
%         xlabel(x_label);
%         colorbar
             %hi = subplot(1,2,1); 

 %            h = pcolor(yy(2:ny-1),xx(2:nx-1),rad_energy(2:nx-1,2:ny-1));
             %set(h, 'EdgeColor', 'none');
             %x_label = sprintf('mu =(%f, %f)',mu(1,1),mu(1,2));
%             time_title = sprintf('t = %f (s)',time);
             %xlabel(x_label);
 %            title(time_title);
  %           colorbar
%                           hi = subplot(1,2,2); 
%              h = pcolor(yy,xx,intensity(:,:,7));
%              set(h, 'EdgeColor', 'none');
%              x_label = sprintf('mu =(%f, %f)',mu(7,1),mu(7,2));
%              time_title = sprintf('t = %f (s)',time);
%              xlabel(x_label);
%              title(time_title);
%              colorbar

            %Plot each angular intensity (recall, ntheta must be even)
               for j=1:na/2
                   hi = subplot(2,3,j); 
                   h = pcolor(xx,yy,intensity(:,:,j)');
                   set(h, 'EdgeColor', 'none');
                   x_label = sprintf('mu =(%f, %f)',mu(j,1),mu(j,2));
                   time_title = sprintf('t = %f (s)',time);
                   xlabel(x_label);
                   title(time_title);
                   colorbar
               end

%Section 5.2 tests
%this is a projection!
%there are boundary effects with Dirichlet bcs. Dont use 2,2 index
% for i=1:na
%     if sqrt(1-mu(i,1)^2 -mu(i,2)^2)> 0.5
%         %mu_z = 0.88
%         quiver(0,0,intensity(nx/2,ny/2,i).*mu(i,1), intensity(nx/2,ny/2,i).*mu(i,2),0,'-b','ShowArrowHead','off');
%     else
%         quiver(0,0,intensity(nx/2,ny/2,i).*mu(i,1), intensity(nx/2,ny/2,i).*mu(i,2),0,'-k','ShowArrowHead','off');
%     end
% hold on;
% end
% hold off; 
%             pause(1.0);

end

