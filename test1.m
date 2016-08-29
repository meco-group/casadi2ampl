import casadi.*

p = SX.sym('p');
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');

v = [x;y;z];
f = x^2 + p*100*z^2;
g = z + (1-x)^2 - y;
nlp = struct('x', v, 'p',p, 'f', f', 'g', g);

% Create IPOPT solver object
solver = nlpsol('solver', 'ipopt', nlp);

% Solve the NLP
data = struct('x0' , [2.5 3.0 0.75],... % solution guess
             'p', 3, ...
             'lbg',    0,...           % lower bound on g
             'ubg',    0);             % upper bound on g
         
 casadi2ampl('results2',nlp,data,[false, false, false])