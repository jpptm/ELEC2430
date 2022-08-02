clc
clear all

%number of samples increased to 70 as per specs
n = 0:1:70; 
%Representation for a 8-point MA filter
x = 2 * sin(((pi*n)/10)-(pi/3));
B = 1/8 * ones(8,1);
output  = filter(B,1,x);
%plot of the input to the filter
stem(n,x)
title('Input');
xlabel('n')
ylabel('x(n)')
%plot of the filter function
stem(B)
title('Moving AVerage Filter');
xlabel('n')
 ylabel('B')
%plot of ouput of the filter
stem(n,output)
title('Output');
xlabel('n')
ylabel('y(n)')