% Generate x and y vectors then plot them
%{x = 0:0.1:30;
y = zeros(length(x));

for num=1:1:length(x)
    y(num) = exp(-0.1 * x(num)) * sin(2 * x(num) / 3);
end

plot(x, y)
xlabel('x');
ylabel('y(x)');

% Part 2
x = -2:1:6;
y = zeros(length(x));

y(3) =  1; y(4) = 2; y(5) = 1; y(6) = 0; y(7) = -1;
stem(x, y, 'filled')
xlabel('x')
ylabel('y')
%}

R = 1; C = 1; T = .2; % Constants, with T as our sampling rate
a = -(1 - T/R/C); b = [0 T/R/C];
y0 = 0; x0  = 1;
n = 1:40;
x = ones(1, length(n));
y1 = recur(a, b, n, x, x0, y0);

plot(n*T, y1, "o");
xlabel("Time (s)")
ylabel("y(t)")


% Reproducing fig 2.16
function y = recur(a, b, n, x, x0, y0)
    N = length(a);
    M = length(b)-1;

    y = [y0 zeros(1,length(n))];
    x = [x0 x];

    a1 = a(length(a):-1:1);    % reverses the elements in a
    b1 = b(length(b):-1:1);

    for i=N+1:N+length(n)
    y(i) = -a1*y(i-N:i-1)' + b1*x(i-N:i-N+M)';
    end

    y = y(N+1:N+length(n));
end

function y = recur(a,b,n,x,x0,y0);
%
% y = recur(a,b,n,x,x0,y0)
%  solves for y[n] from:
% y[n] + a1*y[n-1] + a2*y[n-2]... + an*y[n-N]
%   =  b0*x[n] + b1*x[n-1] + ... + bm*x[n-M] 
%
% a, b, n, x, x0 and y0 are vectors
% a = [a1 a2 ... aN]
% b = [b0 b1 ... bM]
% n contains the time values for which the solution will be computed
% y0 contains the initial conditions for y, in order,
%    i.e., y0 = [y[n0-N], y[n0-N+1], ...,y[n0-1]]
%    where n0 represents the first element of n
% x0 contains the initial conditions on x, in order
%    i.e., x0 = [x[n0-M],...,x[n0-1]]
% the output, y, has length(n)
%
N = length(a);
M = length(b)-1;
if length(y0) ~= N,
  error('Lengths of a and y0 must match')
end
if length(x0) ~= M,
  error('Length of x0 must match length of b-1')
end
y = [y0 zeros(1,length(n))];
x = [x0 x]
a1 = a(length(a):-1:1)             % reverses the elements in a
b1 = b(length(b):-1:1)       
for i=N+1:N+length(n),
  y(i) = -a1*y(i-N:i-1)' + b1*x(i-N:i-N+M)';
end
y = y(N+1:N+length(n))