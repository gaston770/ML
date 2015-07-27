% Excercise 1.10 - Learning From Data - Chapter 1

% This function does exp_n experiments, where in each experiment we flip a coin 10 times to count
% the number of heads and store the proportion as the result of the experiment.
% It returns the outcome of the first experiment, the outcome of a random experiment and the
% outcome of the experiment where the proportion ended up being smallest. The last selection
% amounts to choosing a hypothesis that minimises error with training data.
function [v_1, v_rand, v_min] = runExperiment(exp_n)
  coins_v = zeros(1, exp_n);
  for i=1:exp_n
    % <=0.5 - tails, >0.5 - head
    % flip the coins 10 times and count how many heads
    coins_v(i) = sum(rand(1, 10) > 0.5) / 10;
  end
  v_1 = coins_v(1);
  v_rand = coins_v(floor(rand() * exp_n + 1));
  [unused, iv_min] = min(coins_v);
  v_min = coins_v(iv_min);
end

% Function that estimates P[|v-u|>epsilon], based on the experiment results given as input, and 
% plots it alongside hoeffding's bound for an interval of epsilon between 0 and 100.
function [] = calcP(v_ele) 
  estimation = zeros(100, 1);
  for epsilon = 1:100
    epsilon_ = 0 + 0.01 * epsilon;
    % Cases where the difference between experiment and actual median is greater than epsilon.
    pos = 0;
    for i=1:numel(v_ele)
      pos += abs(v_ele(i)-0.5) > epsilon_;
    end
    estimation(epsilon) = pos/numel(v_ele);
  end
  
  hoeffding = zeros(100, 1);
  for epsilon = 1:100
    epsilon_ = 0 + 0.01 * epsilon;
    hoeffding(epsilon) = 2*e^(-2*(epsilon_**2)*10); 
  end
  plot(estimation);
  hold on;
  plot(hoeffding);
  hold off;
  pause;
end


% Run the experiments.
total_iter = 1000
v_1 = zeros(total_iter, 1);
v_rand = zeros(total_iter, 1);
v_min = zeros(total_iter, 1);
for i=1:total_iter
  [v_1(i), v_rand(i), v_min(i)] = runExperiment(100);
end

% Coordinates of the center of the histogram bars.
histo_coords = [0,1,2,3,4,5,6,7,8,9]*0.10 + 0.05;
figure(1)
hist(v_1, histo_coords);
figure(2)
hist(v_rand, histo_coords);
figure(3)
hist(v_min, histo_coords);
pause;

% Plot the graphs of P and Hoeffding's bound
calcP(v_1);
calcP(v_rand);
% This last graph should show us that selecting the hypothesis after getting the data
% violates the bound.
calcP(v_min);

