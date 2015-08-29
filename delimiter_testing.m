[A,delimiter_out,headerlines_out] = importdata('D:\USA_Census2.txt');
%[B,delimiter_out,headerlines_out] = importdata('D:\USA_Census.txt');

% for i= 1:101
%     
%     A.data(i,1) = i
%     
% end

%%
m = 101

x = 1:m;
high = 100;
%Defining Parameters

for t0 = -high:high
    for t1 = -high:high



% Defining Hypothesis

h = t0 + t1 * x;

% Cost Function

cf=0;

        for i=1:m

            cf = cf + (h(1) - A.data(1,2))^2;

        end

J(t0+high+1,t1+high+1) = (1 / (2 * m)) * cf;

    end
end


%%

t0 = 1:high;
t1 = 1:high;
%plot (J,t0,t1)
minimum = min(J);
minimum = min(minimum);

% Minima Coordinates
for i=1:(2*high+1)
    
    for j=1:(2*high+1)
       
        if (J(i,j) == minimum)
            
            disp('Minima Coordinates are: ')
            disp(i);
            disp(j);
            
        end
    end
    
    
end

surf(J);

