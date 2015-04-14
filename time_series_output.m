function [y_out] = time_series_output(y_out,time_out,rad_energy,rad_flux,rad_pressure,v,vCsquare,nx,ny,C,GasMomentum,GasKE,GasTE,photon_momentum)
%TIME_SERIES_OUTPUT Uncomment plots which are supposed to print out line
%plot values every interval timesteps, including the first one

%Need to constantly update function definition for what variables are needed depending on desired
%output
persistent out_count;

if isempty(out_count)
    out_count = 0;
end
out_count = out_count+1;

%Set to 2 if you want to ignore the initial condition step, t=0
start_index = 1; 

%Jiang14 Figure 2
y_out(out_count,1) = rad_flux(nx/2,ny/2,1);
%advective flux
y_out(out_count,2) = v(nx/2,ny/2,1)*(rad_energy(nx/2,ny/2) + rad_pressure(nx/2,ny/2,1,1))/C;
y_out(out_count,3) = GasMomentum(nx/2,ny/2,1)+photon_momentum(nx/2,ny/2,1);
y_out(out_count,4) = GasMomentum(nx/2,ny/2,1);
y_out(out_count,5) = GasTE(nx/2,ny/2)+GasKE(nx/2,ny/2)+rad_energy(nx/2,ny/2);
y_out(out_count,6) = GasTE(nx/2,ny/2)+GasKE(nx/2,ny/2);
y_out(out_count,7) = rad_pressure(nx/2,ny/2,1,1)/rad_energy(nx/2,ny/2);
y_out(out_count,8) = rad_pressure(nx/2,ny/2,2,2)/rad_energy(nx/2,ny/2);
%relativistic advective flux
y_out(out_count,9) = y_out(out_count,2)/(1+vCsquare(nx/2,ny/2));
figure(1);
%hi = subplot(3,1,1);
plot(time_out(start_index:out_count),y_out(start_index:out_count,1)','-',...
    time_out(start_index:out_count),y_out(start_index:out_count,2)','-',...
    time_out(start_index:out_count),y_out(start_index:out_count,9)','-');
legend('F_{r,x}','F_{v,x} O(v/c)','F_{v,x} SR');
xlabel('time (s)');
%hi = subplot(3,1,2);
%plot(time_out(start_index:out_count),y_out(start_index:out_count,3)','--',time_out(start_index:out_count),y_out(start_index:out_count,4)','-');
%xlabel('time (s)');
%legend('p_{x,tot}','p_{x,gas}');

%hi = subplot(3,1,3);
%plot(time_out(start_index:out_count),y_out(start_index:out_count,5)','--',time_out(start_index:out_count),y_out(start_index:out_count,6)','-');
%axis([0 1.0 4.5 7.5]);
%xlabel('time (s)');
%legend('E_{tot}','E_{gas}');
% 
% figure(2);
% axis([0 1.0 0.3 0.4]);
% plot(time_out(start_index:out_count),y_out(start_index:out_count,7)','-k',time_out(start_index:out_count),y_out(start_index:out_count,8)','--k');
% xlabel('time (s)');
% legend('f_{x,x}','f_{y,y}');




%pause(1.0);


     %Figure 1
%     y_out(out_count) = rad_energy(nx/2,ny/2);
%     y_out2(out_count) = temp(nx/2,ny/2)^4;
% % %     figure(1);
%       legend('E_r','T^4');

%      xlabel('time (s)');
%      ylabel('E_r,T^4');
%      title('\sigma_a = 1');


end

