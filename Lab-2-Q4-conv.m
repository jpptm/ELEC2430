%---------Option A----------%

n = 0:1:40;
hn = exp(-0.01.*n).*sin(0.01.*n);
isn = 1 - cos(0.01 .* n);
Io = conv(hn, isn);
n = 0:1:80;

stem(n, Io);
title 'Discrete Time Approximation'
xlabel n
ylabel Io[n]


%---------Option B----------%
n = 0:1:40;
hn = exp(-0.01.*n).*sin(0.01.*n);
isn = 1 + sin(0.01 .* n);
Io = conv(hn, isn);
n = 0:1:80;

stem(n, Io);
title 'Discrete Time Approximation'
xlabel n
ylabel Io[n]
