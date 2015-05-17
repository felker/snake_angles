figure(12);
j=6;
l=1;
abs_error = abs(intensity_analytic(is:ie,js:je,j,l) -intensity(is:ie,js:je,j,l));
combination = abs(intensity_analytic(is:ie,js:je,j,l) + intensity(is:ie,js:je,j,l));
h= pcolor(xx,yy,abs_error');
set(h, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title('| I - I_s |');

%collect statistics about the solution error
% max,mean,rms erros
percents = zeros(size(abs_error,1),size(abs_error,2));
err_stats = [max(abs_error(:)),mean(abs_error(:)),rms(abs_error(:))]
for i=1:size(abs_error,1)
    for j=1:size(abs_error,2)
    %cant divide by zero if both analytic and numeric solution are zero!
    
    %we are getting very small numerical solutions that mess up the error
    %when the analytic solution is 0.0
        if intensity_analytic(i+num_ghost,j+num_ghost,6,1) == 0
            percents(i,j) = 0.0;
        elseif combination(i,j) < 1 %now, i choose to ignore points not along beam
            percents(i,j) = 0.0;
        else
            percents(i,j) = 200*abs_error(i,j)/combination(i,j);
        end
    end
end

err_percents = [max(percents(:)),mean(percents(:)),rms(percents(:))]
figure(13);
h= pcolor(xx,yy,percents');
set(h, 'EdgeColor', 'none');
colorbar;
title('Percent error inside beam')
xlabel('x');
ylabel('y');