%---------Option A----------%
% A
T = 0.01; % Assign sampling period
n = 0:T:20; % Define sampling vector

% Evaluate functions
hn = exp(-n).*sin(n);
isn = 1 - cos(n);

% Take convolution then truncate to match length
Io = T.*conv(hn, isn);
Io = Io(1:length(n));

% Plot
figure
stem(n, Io);
grid
title 'Discrete Time Approximation'
xlabel n
ylabel Io[n]


%---------Option B----------%
T = 0.01; % Assign sampling period
n = 0:T:20; % Define sampling vector

% Evaluate functions
hn = exp(-n).*sin(n);
isn = 1 + sin(n);

% Take convolution then truncate to match length
Io = T.*conv(hn, isn);
Io = Io(1:length(n));

% Plot
figure
stem(n, Io);
grid
title 'Discrete Time Approximation'
xlabel n
ylabel Io[n]