%% Basic Perceptron Learning Algorithm
% Implementation of Learning From Data book for the Perceptron learning algorithm.
% It also graphs the decision boundary of the perceptron on each itearation, so you can see
% what 'learning' looks like for it.

% This is the input data matrix, or feature matrix. Each row corresponds to a training
% sample, and the columns of a row represent the features that sample. 
X = [[1 1],
     [1 4],
     [2 5],
     [1 0],
     [2 1],
     [2 2]];

% This vector represents the class of each sample data. A perceptron can only separate the state
% space in two, so we have two classes to choose from: -1 and +1. The first three samples belong
% to class +1 and the last three samples belong to class -1.
y = [1 1 1 -1 -1 -1];

m = rows(X);

% For plotting the samples we'll divide them in their respective groups.
A = X(1:3, :);
B = X(4:6, :);

% This function evaluates a sample with our hypothesis, returning -1 or +1 in 'val'. Since
% samples don't have bias column in them, needed to allow a hypothesis to learn a border that
% doesn't intersect the (0, 0) coordinates, we add it to the sample. 
function [val] = evalTheta(X, theta)
  X = [1 X];
  val = sign(theta * X');
end

% This function is used to plot the border of the hypothesis. That is, where is the function value
% equal to 0. To one side of it, we'll have the class +1 and to the other side the class -1.
% As any function plotted in 2D grid, it takes only one coordinate, the x axis position.
function [val] = getBorder(x, theta)
  % From sign(theta(1) + theta(2) * x_1 + theta(3) * x_2) we clear x_2 (the y coordinate in our 
  % samples plot) and come up with this formula.
	val = -(theta(1) + theta(2) * x) / theta(3);
end

% Our learning algorithm starts with a random hypothesis and has a rate of adjustment set to 0.01.
% You should try changing that and seeing how it affects the learning speed and hypothesis 
% complexity. A very low value will give use very low convergence speed, but it will also navigate
% more complex hypothesies. In escence what we are doing is trying to find a linear combination of
% the hypothesis and the samples such that its resulting vector hypothesis classifies the samples
% correctly. 

% The algorithm all it does is to navigate the space of possible solutions trying to fit,
% at each step of the algorithm, an example that is currently not being classified correctly.

% If we consider the iterations as t = 1, 2, 3, 4... the key idea is the following: 
%  % Given a misclassified sample (x, y) update w like this "w(t+1)=w(t)+x(t)*y(t)"
% We can see that such sample x(t) has the following properties:
%   * w(t)'*x(t)*y(t) < 0 (i.e. it is missclassified) -> this follows from h(x(t)) != y(t)
%   * w(t+1)'*x(t)*y(t) > w(t)'*x(t)*y(y) (i.e. closer to 0) -> this follows from the update rule

done = false;
theta = rand(1, 3);
rate = 0.01;
while !done
  updated = false;
  for i=1:m
    % Find only the first unclassified sample and update weights. We could do it for all 
    % misclassified samples too.
    if sign(evalTheta(X(i, :), theta)) != y(i) and (!updated);
      % Bias term is never updated.
      theta = theta + rate * [0 X(i, :)] * y(i);
      updated = true;
    end
  end
  
  % Plot the data set and the learning lines
  plot(A(:, 1), A(:, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);
  hold on;
  plot(B(:, 1), B(:, 2), 'bx', 'MarkerSize', 10, 'LineWidth', 1.5);
  hold off;
  hold on;
  plot ([0, 3], [getBorder(0, theta), getBorder(3, theta)], '--', 'MarkerSize', 10, 'LineWidth', 1.5 );
  hold off;
  sleep(0.3);
  if !updated
    done = true;
  end
end
fprintf("Final Perceptron's weights\n");
theta
pause;
