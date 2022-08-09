clc
clear all

%-----------------OPTION A-----------------%
% Index
n = 0:50; 
% Unit step signal
u = (n >= 0); 
% h[n] 
h = 2*0.6.^n.*u; 

% A
x1 = u; 
% Output to Response A
a_y = conv(h,x1); 

% B
x2 = cos(n*pi/4).*u; 
% Output to Response B
b_y = conv(h,x2);

% C
x3 = u + cos(n*pi/4).*u; 
% Output To Response C
c_y = conv(h,x3); 

%------------Response Plot-------------%
% Plot of response A
stem(a_y)
% Plot of response B
stem(b_y) 
% Plot of response C
stem(c_y) 
% Plot of response a+b=c satisfies linearity check
stem(a_y + b_y)

Result = 'From Observing the Plots you can realise that the response of c_y is the same as response of  a_y + b_y hence system is linear';
disp(Result)


%-----------------OPTION B-----------------%
% Index
n = 0:50; 
% Unit step signal
u = (n >= 0); 
% h[n] 
h = 0.5.^n + 1.*u; 

% A
x1 = u; 
% Output to Response A
a_y = conv(h,x1); 

% B
x2 = cos(n*pi/4).*u; 
% Output to Response B
b_y = conv(h,x2);

% C
x3 = u + cos(n*pi/4).*u; 
% Output To Response C
c_y = conv(h,x3); 

%------------Response Plot-------------%
% Plot of response A
stem(a_y)
% Plot of response B
stem(b_y) 
% Plot of response C
stem(c_y) 
% Plot of response a+b=c satisfies linearity check
stem(a_y + b_y)

Result = 'From Observing the Plots you can realise that the response of c_y is the same as response of  a_y + b_y hence system is linear';
disp(Result)
