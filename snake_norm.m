function [ norm ] = snake_norm(v,alpha,beta)
%snake_norm Compute the norm of column vector
%   with the snake metric tensor given alpha and beta
    dims = size(v);
    assert( dims(2) == 1);
    switch dims(1)
        case 2 %induced metric on spacelike 2D hyperplane
             norm = (alpha*v(1)).^2 - 2*beta*v(1)*v(2) + v(2)^2 ;
        case 3 %induced metric on spacelike 3D hyperplane
            norm = (alpha*v(1)).^2 - 2*beta*v(1)*v(2) + v(2)^2 + v(3)^2;
        case 4 %full metric 
            norm = -(v(1)^2) + (alpha*v(2)).^2 - 2*beta*v(2)*v(3) + v(3)^2 + v(4)^2;
        otherwise
            error('Vector has incorrect dimensions!');
    end
end

%NEED TO FIX INDEXING WHEN DEALING WITH 4D vectors

