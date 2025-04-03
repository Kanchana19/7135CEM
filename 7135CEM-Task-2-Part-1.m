% Create a Fuzzy Logic Controller (FLC) for Fan Speed Control

% Define fuzzy logic system
fis = mamfis('Name', 'FanSpeedControl');

% Define input variable: Temperature
fis = addInput(fis, [0 40], 'Name', 'Temperature');
fis = addMF(fis, 'Temperature', 'trimf', [0 0 20], 'Name', 'Cold');
fis = addMF(fis, 'Temperature', 'trimf', [10 20 30], 'Name', 'Warm');
fis = addMF(fis, 'Temperature', 'trimf', [20 40 40], 'Name', 'Hot');

% Define input variable: Humidity
fis = addInput(fis, [0 100], 'Name', 'Humidity');
fis = addMF(fis, 'Humidity', 'trimf', [0 0 50], 'Name', 'Low');
fis = addMF(fis, 'Humidity', 'trimf', [25 50 75], 'Name', 'Medium');
fis = addMF(fis, 'Humidity', 'trimf', [50 100 100], 'Name', 'High');

% Define output variable: Fan Speed
fis = addOutput(fis, [0 100], 'Name', 'FanSpeed');
fis = addMF(fis, 'FanSpeed', 'trimf', [0 0 50], 'Name', 'Slow');
fis = addMF(fis, 'FanSpeed', 'trimf', [25 50 75], 'Name', 'Medium');
fis = addMF(fis, 'FanSpeed', 'trimf', [50 100 100], 'Name', 'Fast');

% Define fuzzy rules
ruleList = [
    1 1 1;
    1 2 2;
    1 3 2;
    2 1 1;
    2 2 2;
    2 3 3;
    3 1 2;
    3 2 3;
    3 3 3
];

fis = addRule(fis, ruleList);

% Show FIS structure
disp(fis);

% Plot membership functions
figure;
subplot(3,1,1); plotmf(fis, 'input', 1);
title('Temperature Membership Functions');
subplot(3,1,2); plotmf(fis, 'input', 2);
title('Humidity Membership Functions');
subplot(3,1,3); plotmf(fis, 'output', 1);
title('Fan Speed Membership Functions');

% Evaluate FIS
input_values = [30 60]; % Example input: Temperature = 30, Humidity = 60
output = evalfis(fis, input_values);

% Display result
fprintf('Fan Speed Output: %.2f\n', output);