% Function to calculate the decoding matrix and efficiency of an encoding 
% scheme given the encoding matrix.  Also returned are the rank and
% condition of the matrix (R and C respectively).
%
% function [A_p, E, R, C] = AssessEncoding(A);

function [A_p, E, R, C] = AssessEncoding(A);

A_p = pinv(A);

E = 1./sqrt(sum(A_p.^2,2)*size(A_p,2));

R = rank(A);

C = cond(A);


disp(['Matrix rank is: ' num2str(R)])
disp(['Matrix condition is: ' num2str(C)])
disp('SNR efficiency for each component is: ')
disp(E)