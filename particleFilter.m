% John Lawler
% Lab 6 - Particle Filter
clc;
clear all;
close all;
%%
%% Read in data
file = fopen('magnets-data.txt', 'r');
format = '%f %f %f';
Data = textscan(file, '%f %f %f');
fclose(file);
data_mat = cell2mat(Data);
actual_pos = transpose(data_mat(:, 1));
actual_vel = transpose(data_mat(:, 2));
sensor_val = transpose(data_mat(:, 3));
data_size = size(sensor_val)
count = data_size(2);
time = 1:size(data_mat);

%% Initialize variables

% M is # of particles chosen to represent the unkown distribution p(x|y)
% Increase M as the dimsionality and shape complexity of p(x|y) increases
M = 100;
init_weight = 1/M;
init_state = [1; 0];
particles = {init_state, init_weight};
updated_weights = [];

% Initialize particles and weights to a cell array of M rows (# particles)
for i=2:M
    particles(end+1, :) = {init_state, init_weight};
end

% how to retrieve a state
% state1 = cell2mat(particles(:, 1))
% state1(1)

% how to retrieve a weight
% weight1 = cell2mat(particles(1, 2))


sigma_a = 2^(-4);
sigma_m = 4.0;
sigma_n = 2^(-8);


get_plot_points = 15;



%% Predict-Update Cycle

% Loop for amonut of sensor readings
for i=1:count 
  
  % Record previous points for plotting
  if i == get_plot_points
     plot_pts_tm1 = particles(:, 1);
     plot_wts_tm1 = particles(:, 2);
  end
  
  % Step 1: Propogate particles through state transitions
  for j=1:M
      % Retrieve a particle's previous state variables
      x_tm1 = cell2mat(particles(j, 1));
      
      % Predict guesses based on gaussian noise
      Xt_m = state_eq(x_tm1(1), x_tm1(2), N(0, sigma_a));
      
      % Set new state variables back to particle list
      particles(j, 1) = {Xt_m};
  end
  
  % Record points for plotting
  if i == get_plot_points
     plot_pts_t = particles(:, 1);
     plot_wts_t = particles(:, 2);
  end
  
  
  % Step 2: Retrieve Yt
  yt = sensor_val(i);
  
  % Step 3: Update the weight for each particle using the measurements
  for j=1:M
      % Retrieve previous weight value w_(t-1)
      w_tm1 = cell2mat(particles(j, 2));
      
      % Retrieve xt_m
      xt_state = cell2mat(particles(j, 1));
      xt_m = xt_state(1);
      
      % Calculate p(yt | xt(m) )
      yt_m = ideal_meas(xt_m, sigma_m);
      p_dis = distr(yt, yt_m, sigma_n);
      
      % Calculate ~w
      updated_weights(j) = w_tm1 * p_dis;
  end 
  
  % Step 4: Normalize the weights
  weight_sum = sum(updated_weights);
  % updated_weights(1)/weight_sum
  for j=1:M
      particles(j, 2) = {updated_weights(j)/weight_sum};
  end
  
  % Step 5: compute desired output (expected value)
  for j=1:M
      % Retrieve xt_m and wt_m
      xt_state = cell2mat(particles(j, 1));
      xt_m = xt_state(1);
      xdt_m = xt_state(2);
      wt_m = cell2mat(particles(j, 2));
      
      E(j) = xt_m * wt_m;
      Ed(j) = xdt_m * wt_m;
  end
  
  
  exp_pos = sum(E);
  output_pos(i) = exp_pos;
  
  exp_vel = sum(Ed);
  output_vel(i) = exp_vel;
    
  % Step 6: Find ESS and check for resample
  weights = cell2mat(particles(:,2));
  cv = calculate_CV(M, weights);
  ESS = M/(1+cv);
  
  % Resample
  if ESS < 0.5*M 
      pstates = particles(:,1);
      Q = cumsum(weights);
      t=rand(M+1);
      T=sort(t);
      T(M+1)=1.0;
      k=1;
      r=1;
      while (k <= M)
          if T(k) < Q(r)
              Index(k)=r;
              k=k+1;
          else
              r=r+1;
          end
      end
      for c=1:M
         newP = cell2mat(pstates(Index(c)));
         newW = 1/M;
         
         particles(c, :) = {newP, newW};
      end
 
      % Record points for plotting
      if i == get_plot_points
         pts_r = particles(:, 1);
         wts_r = particles(:, 2);
      end 
  end
 
end


%% Data Plot of Actual and Estimate Position full range
figure
hold on
grid off
axis([0 1100 -20 20]);
%set(gca, 'FontSize', 14)
a = plot(time, actual_pos, 'x');
y = plot(time, output_pos);
set(y, 'Color', '#f53b02');
set(a, 'Color', '#BAB6B3');
xlabel('Time (s)')
ylabel('Position')
leg = legend('Actual','Estimate','Location', 'southeast');
set(leg, 'FontName', 'Helvetica');

%% Data Plot of Actual and Estimate Velocity full range
figure
hold on
grid off
axis([0 1100 -3 3]);
%set(gca, 'FontSize', 14)
a = plot(time, actual_vel, 'x');
y = plot(time, output_vel);
set(y, 'Color', '#f53b02');
set(a, 'Color', '#BAB6B3');
xlabel('Time (s)')
ylabel('Position')
leg = legend('Actual','Estimate','Location', 'southeast');
set(leg, 'FontName', 'Helvetica');


%% Data Plot of Actual and Estimate Position for time 0 to 100
figure
hold on
grid off
axis([0 100 -1.5 2.5]);
%set(gca, 'FontSize', 14)
a = plot(time, actual_pos, 'x');
y = plot(time, output_pos);
set(y, 'Color', '#f53b02');
set(a, 'Color', '#BAB6B3');
xlabel('Time (s)')
ylabel('Position')
leg = legend('Actual','Estimate','Location', 'southeast');
set(leg, 'FontName', 'Helvetica');

%% Data Plot of a set of particles for a given time, choosing time = 50, weights vs position
figure
hold on
grid off
%axis([0 100 -1.5 2.5]);
%set(gca, 'FontSize', 14)
wts = cell2mat(plot_wts_t);
for i=1:M
    xt_state = cell2mat(plot_pts_t(i));
    xt_m = xt_state(1);
    pts(i) = xt_m;
end

a = plot(wts, pts, 'x');
%set(y, 'Color', '#f53b02');
set(a, 'Color', '#2D2B2A');
xlabel('Weight')
ylabel('Position')
leg = legend('Particles','Estimate','Location', 'southeast');
set(leg, 'FontName', 'Helvetica');


%% Data Plot of Set of particles before and after time t = 15
figure
hold on
grid off
%axis([0 100 -1.5 2.5]);
%set(gca, 'FontSize', 14)
wts_tm1 = cell2mat(plot_wts_tm1);
for i=1:M
    xt_state = cell2mat(plot_pts_tm1(i));
    pts_tm1(i) = xt_state(1);
end

b = plot(wts_tm1, pts_tm1, 'x');
a = plot(wts, pts, 'o');
set(b, 'Color', '#f53b02');
set(a, 'Color', '#2D2B2A');
xlabel('Weight')
ylabel('Position')
leg = legend('Particles at t','Particles at t-1','Location', 'southeast');
set(leg, 'FontName', 'Helvetica');

%% Resample plot of time t = 15
figure
hold on
grid off
%axis([0 100 -1.5 2.5]);
%set(gca, 'FontSize', 14)
r_wts = cell2mat(wts_r);
for i=1:M
    xt_state = cell2mat(plot_pts_t(i));
    pts(i) = xt_state(1);
    
    xt_state = cell2mat(pts_r(i));
    r_pts(i) = xt_state(1);
end

a = plot(wts, pts, 'x');
b = plot(r_wts, r_pts, 'o');

set(b, 'Color', '#f53b02');
set(a, 'Color', '#2D2B2A');
xlabel('Weight')
ylabel('Position')
leg = legend('Before','Afer','Location', 'southeast');
set(leg, 'FontName', 'Helvetica');


%% Measurement Plot
figure
grid off
axis([0 1100 -20 20]);
d = plot(time, sensor_val, '.');
xlabel('Time (s)')
ylabel('Measured Position')

% Generate random noise based on gaussian distribution
function at = N(a, b)
    at = normrnd(a,b);
end

% Calculate state transitions
function Xt = state_eq(xt, xd, at)
    T = 1;
    % x_(t+1)
    xtp1 = @(xt, xd) (xt + xd*T);
   
    % xdot_(t+1)
    if xt < -20 
        xdp1 = 2;
    elseif (-20 <= xt) && (xt < 0)
        xdp1 = xd + abs(at);
    elseif (0 <= xt) && (xt <= 20)
        xdp1 = xd - abs(at);
    elseif xt > 20
        xdp1 = -2;
    end
    
    Xt = [xtp1(xt, xd); xdp1];
end

% Calculate ideal measurement of a given particle
function yt_m = ideal_meas(xt_m, sig_m)
    xm1 = -10;
    xm2 = 10;
    temp = (1/(sqrt(2*pi)*sig_m))*exp(-((xt_m - xm1)^2)/(2*(sig_m^2)));
    yt_m = temp + (1/(sqrt(2*pi)*sig_m))*exp(-((xt_m - xm2)^2)/(2*sig_m^2));
end

% Calculate probability distribution for a given particle
function p = distr(yt, yt_m, sig_n)
    p = (1/(sqrt(2*pi)*sig_n))*exp(-((yt_m - yt)^2)/(2*sig_n^2));
end

% Calculate coefficient of variation for resampling check
function cv = calculate_CV(M, weights)
   for i=1:M
       val(i) = (M*weights(i) - 1)^2;
   end
   cv = (1/M)*sum(val);
end



