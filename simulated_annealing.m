clear;clc
%% Ackley_2D
Ackley_2D = @(x1,x2,a,b,c) -a.*exp(-b.*sqrt(x1.^2 + x2.^2))-exp((cos(c.*x1)+cos(c.*x2))./2) + a + exp(1);
frist = -40; last = 40; 
x1 = frist:0.1:last ; x2 = frist:0.1:last;

a = 10 ;
b = 0.2;
c= pi/2;
%plot(x1,x2,Ackley_2D(x1,x2,a,b,c))

[X,Y] = meshgrid(x1,x2);
Z = Ackley_2D(X,Y,a,b,c);
figure(1);
colorbar;
surf(X,Y,Z,"LineStyle","none");
title('2D Ackley function')
xlabel('x_{1}')
ylabel('x_{2}')
%axis equal;
%view(0,90);

%% Paraboloid
para = @(x1,x2) x1.^2 + x2.^2 ;
P = para(X,Y);
figure(2);
colorbar;
surf(X,Y,P,"LineStyle","none");
title('Paraboloid')
xlabel('x_{1}')
ylabel('x_{2}')
axis equal;
view(0,90);
%% simulated annealing
bounds = [frist last;frist last];
%define the total iterations
n_iterations = 1000;
step_size = 0.1;
%initial temperature
temp = 10;
%%%%
m = ones(n_iterations,1);
best = bounds(:,1) + rand(length(bounds),1) .* (bounds(:, 2) - bounds(:, 1));
%best = [1 1];
%best_eval = Ackley_2D(best(1),best(2),a,b,c);
best_eval = para(best(1),best(2));
curr = best; curr_eval = best_eval;
for i = 1:1:n_iterations
    candidate = curr + randn(length(bounds),1) .* step_size;
    %candidate_eval = Ackley_2D(candidate(1),candidate(2),a,b,c);
    candidate_eval = para(candidate(1),candidate(2));
    if candidate_eval < best_eval
        best  = candidate ;
        best_eval = candidate_eval;
        disp(best_eval);
    end
    %Energy change
    diffE = candidate_eval - curr_eval;
    %temperature
    t = temp / (i + 1);
    %t = temp;

    metropolis = exp(-diffE / t);
    if diffE < 0 || rand() < metropolis
        curr  = candidate; 
        curr_eval = candidate_eval;
        disp(curr_eval);
        disp(i);
    end
    %disp(i);
    m(i,1) = i;
    m(i,2) = best_eval;
end
%disp(i);
title_graph = 'Evaluation of Paraboloid';
%title_graph = '2D Ackley function';
X = sprintf('x1 = %d , x2 = %d and z = %d.',best(1) ,best(2),best_eval);
disp(X)

figure(3);
plot(m(:,1),m(:,2));
title(title_graph);
xlabel('iterations')
ylabel('Evaluation f(x)')

