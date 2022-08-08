clc
clear all

% 8 point MA filter
h = (ones(1,8))./8;

% Overshooting of n to 70 samples instead of 40
n = 0:1:70;

% input
x = 2 * (sin((pi*n / 10) - (pi / 3)));

% Create a matrix and fill it with zeros
N1 = length(x);
N2 = length(h);
x = [x,zeros(1,N2)];
h = [h,zeros(1,N1)];

% Moving Average filter output
for i = 1: N1 + N2 - 1
    y(i) = 0;
    for j = 1:N1
        if((i-j+1)>0)
        y(i)=y(i)+(x(j)*h(i-j+1));
        end
    end
end

n = 0:1:length(y)-1;

% Plot of Filter input
stem(x);
grid;
xlabel('n');
ylabel('y(n)')
title('Filter Input')

% Plot of Output Response
stem(n,y);
grid;
xlabel('n');
ylabel('y(n)')
title('Option A-8P-MAF Output')



%-------------------OPTION B-----------------%
h = (ones(1,8))./8;

% Overshooting of n to 90 samples instead of 60
n = 0:1:90;

% input
x = 2 * (cos((pi*n / 10) - (pi / 3)));

% Create a matrix and fill it with zeros
N1 = length(x);
N2 = length(h);
x = [x,zeros(1,N2)];
h = [h,zeros(1,N1)];

% Moving Average filter output
for i = 1: N1 + N2 - 1
    y(i)=0;
    for j = 1:N1
        if((i-j+1)>0)
        y(i)=y(i)+(x(j)*h(i-j+1));
        end
    end
end

n=0:1:length(y)-1;

% Plot of Filter input
stem(x);
grid;
xlabel('n');
ylabel('y(n)')
title('Filter Input')

% Plot of Output Response
stem(n,y);
grid;
xlabel('n');
ylabel('y(n)')
title('Option B-8P-MAF Output')
