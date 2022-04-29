function [results] = Create_Stream_Function(results)
%CREATE_STREAM_FUNCTION Creates Stream function from the velocity field

% Define u and v component of the velocity field
u  = results.u;
v  = results.v;
xi = results.xi; 

psi = poisolv()

end

