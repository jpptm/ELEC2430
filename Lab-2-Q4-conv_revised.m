%---------Part A---------%
clear;

plotTo = 20;
syms t tau;

h = exp(-t)*sin(t);
IsA = (1-cos(t))*(heaviside(t)-heaviside(t-plotTo));
IsB = (1+sin(t))*(heaviside(t)-heaviside(t-plotTo)); 

% the calculation for the terms inside each integral is performed here
toIntegrateA = subs(h,tau)*subs(IsA,t-tau);
toIntegrateB = subs(h,tau)*subs(IsB,t-tau);

% the integrations for each convolution are calculated and printed to display
IoA = int(toIntegrateA,tau,0,t)
IoB = int(toIntegrateB,tau,0,t)

% plot for option A
fplot(t, IoA, [0 plotTo]); 
title('Continuous A: Io(t) = Is(t) + h(t)');
ylabel('Io(t)');
xlabel('t')
ylim([0 1])

% plot for option B
fplot(t, IoB, [0 plotTo]);
title('Continuous B: Io(t) = Is(t) + h(t)');
ylabel('Io(t)');
xlabel('t')
ylim([0 1])

%-----------Part B----------%

plotTo = 20; % value the convolution will occur to 
T = 0.01; % sample rate
n = 0:T:plotTo; % discrete values to convolute over

% definitions for h and the Is values of part A and B
h = exp(-n).*sin(n);
IsA = (1-cos(n));
IsB = (1+sin(n));

% performing convolutions and scaling for parts A and B
IoA = conv(h, IsA).*T;
IoB = conv(h, IsB).*T;

% chopping of convolutions to correct dimension
IoA = IoA(1:length(n));
IoB = IoB(1:length(n));

% displaying the graphs for the discrete results
stem(n, IoA, 'linestyle','none'); 
ylim([0 1]); 
title('Discrete A: Io[n] = Is[n] + h[n]');
ylabel('Io[n]');
xlabel('n');

stem(n, IoB, 'linestyle','none'); 
ylim([0 1]); 
title('Discrete B: Io[n] = Is[n] + h[n]');
ylabel('Io[n]');
xlabel('n');

%---------Part C---------%

Looking at the plots for both the continuous and discrete-time methods reveals that these plots look identical to
each other when viewed at this scale. This is because the small sample rate of 0.01 for the discrete convolution
forms an excellent approximation of the solution for the continuous version.