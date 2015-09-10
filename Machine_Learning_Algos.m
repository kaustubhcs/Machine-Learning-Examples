%% Machine Learning Algorithm
% This matlab program helps us understand the use of Hypothesis and cost
% functions on any kind of data _Currently Implemented for USA Census Data_
% and thus applies machine learning algorithms.

%% Clearing Out Screen
% This part of the code snippet should be deleted if using as a 
% function
clear all
clc


%% Load Data


[A,delimiter_out,headerlines_out] = importdata('D:\USA_Census2.txt');

%[B,delimiter_out,headerlines_out] = importdata('D:\USA_Census.txt');


%% Defining Parameters
m = 100;

x = 1:m+1;
high = 100;
data_type=2;


%% Hypothesis and Cost Function Generation

for t0 = -high:high
    for t1 = -high:high



%% Defining Hypothesis

h = t0 + t1 * x;

%% Cost Function

cf=0;

        for i=1:m

            cf = cf + (h(1,m) - A.data(m,data_type))^2;

        end

J(t0+high+1,t1+high+1) = (1 / (2 * m)) * cf;

    end
end


%% Defining Sub parameters

t0 = 1:high;
t1 = 1:high;
%plot (J,t0,t1)
minimum = min(J);
minimum = min(minimum);

%% Minima Coordinates
k=1;
for i=1:(2*high+1)
    
    for j=1:(2*high+1)
       
        if (J(i,j) == minimum)
            
            minimas(k,1) = i - high - 1;
            minimas(k,2) = j - high - 1;
            
            k = k + 1;
            
        end
    end
    
    
end

%%
% *%% Plotting*

axes = 1; 

figure(1);

surf(J);
title('Cost Function');
figure(2);

subplot(2,2,1);
plot(A.data(:,data_type));
title('Data Plot');

subplot(2,2,2);
[msize_1,msize_2] = size(minimas);

for i=1:msize_1
    
    for j = 1:m+1
        
        min_ht(i,j) = minimas(i,1) + minimas(i,2) * x(i);

    
    
    end
    
end


aux_ht = 4 + 0.000004 * x;
plot(x,min_ht(1,:),x,min_ht(2,:));
title('Hypothesis Plot');

subplot(2,2,3);

% 
%%
% Scaling

plot(x,A.data(:,data_type),x,min_ht(1,:),x,min_ht(2,:));
title('Data + Hypothesis Plot');

subplot(2,2,4);
plot(minimas);
title('Minima Plot');

