%This set of functions fit IP model to a experimental data by
%optimization.
%Initial guesses of the parameters are taken from experimental data (for 
%initial conditions).
%
%@Chanwoo Chun, <cc2465@cornell.edu>

function data_output = fitIP(data,st)

%convert yexp (mm) to (cm) by dividing by a factor of 10. This change is not going to be saved.
data.com = data.com;
data.vel3D = data.vel3D;

%Comment out below two lines of code to enable fitting to a pure tripod
%stance phases. Uncomment if you want to fit to entire tripod stance phase.
data.pureTriStarts = 1;
data.pureTriEnds = size(data.com,1);

% Construct fixed structure
fixed.g      = 9807; %mm/s^2
fixed.M      = data.source.weight/1000; %gram
expLegLength = data.source.legLength;

%Definitions and units.
%gammaA = ka/(mgR)
%Ka unit  = kg(m/s)^2 = N*m*rad^-1 = 10^3gram*(10^3mm/s)^2 = 10^9gram*(mm/s)^2
% kg(m/s)^2 / 10^9 gram(mm/s)^2 = 1
%gammaA unit = gram(mm/s)^2/(gram * mm/s^2 *mm) = 1
%Omega unit = rad/s 

%Set up unconstrained minimization
objective = @(X) objectiveFunc(X, data, fixed, st);

A = []; b = []; % no linear inequality constraints
Aeq = []; beq = []; % no linear equality constraints

pureTriStarts = data.pureTriStarts(1);
xzPosInit = data.com(pureTriStarts,:);
xzVelInit = data.vel3D(pureTriStarts,:); 

AP = median(data.com(:,1));% Initial anchor point
travelLength = data.com(end,1)-data.com(1,1);
APlb = data.com(1,1)-travelLength/3; %define lower bound
APub = data.com(end,1)+travelLength/3; %define upper bound
R0 = sqrt(xzPosInit(2)^2+(AP-xzPosInit(1))^2);% Initial leg length
%Calculate dRdt0 and omega
A = -AP+xzPosInit(1);
B = xzPosInit(2);
legUnitVector = [A, B]/norm([A, B]);
dRdt0 = sum(legUnitVector.*xzVelInit); %dot product
tanVel0 = xzVelInit-dRdt0*legUnitVector;
tanSpeed0 = sqrt(tanVel0(1)^2+tanVel0(2)^2);
omega = tanSpeed0/R0;

%Initial guess of the parameter values
params0 = [R0; AP; omega]; %R0 AP omega
lb = [R0*0.1; APlb; omega*0.25]; % lower bound for parameters
ub = [R0*5; APub; omega*5]; % upper bound for parameters Hooke's
lbTemp = min([lb ub],[],2);
ubTemp = max([lb ub],[],2);
lb = lbTemp;
ub = ubTemp;

% Option 1:
%options = optimoptions('fmincon','Algorithm','sqp','Display','iter', ...
%   'GradObj','off','MaxIterations',1500,'MaxFunctionEvaluations',3000);
%
% Option 2:
options = optimoptions('fmincon','Algorithm','interior-point','Display','iter', ...
     'GradObj','off');
problem = createOptimProblem('fmincon','x0',params0,'objective',objective,...
    'lb',lb,'ub',ub,'options',options);

gs = GlobalSearch('NumTrialPoints', 2000);
%try
    [reconstructed, f] = run(gs, problem);
% catch
%     %If ran in to an error during optimization, save NaN to the parameters.
%     data_output=data;
%     return
% end

%Assign optimized parameter values.
R      = reconstructed(1);
AP      = reconstructed(2);
omega   = reconstructed(3);
APpct = reconstructed(2)/travelLength*100;

%Solve ODE to save optimized solutions.
thetasim_interp = getCoordinates(reconstructed, data, fixed, st);
% Compute predicted trajectory in x and z directions
ysim(:,1) = R.*(sin(thetasim_interp(:,1))-sin(thetasim_interp(1,1)))+data.com(data.pureTriStarts(1),1);
ysim(:,2) = R.*(cos(thetasim_interp(:,1)));
% Raw solution (angle, angular speed, leg length, and leg speed)
rawSolution = thetasim_interp;

data.IP.R = R;
data.IP.AP = AP;
data.IP.APpct = APpct;
data.IP.omega = omega;
data.IP.f = f;
data.IP.ysim = ysim;
data.IP.rawSolution = rawSolution;
data.IP.movesBack = false;
if any(ysim(:,2)<=0)
        data.IP.movesBack = true;
end
data_output=data;
end

%%

function [ f, g ] = objectiveFunc(params, data, fixed, st)
%Compute the score
%   params: [Rnat R0 Ks AP dRdt0 omega]
%   yexp: an Ntx2 vector of data in x and z directions
%   fixed: a structure of fixed quantities
pureTriStarts = data.pureTriStarts(1);
pureTriEnds = data.pureTriEnds(1);

% Load data
yexp = data.com;

% Compute the score of current params
% First, solve for theta and R using ODE
thetasim_interp = getCoordinates(params, data, fixed, st);
R = params(1);
% Compute predicted trajectory in x and z directions
ysim(:,1) = R.*(sin(thetasim_interp(:,1))-sin(thetasim_interp(1,1)))+yexp(pureTriStarts,1);
ysim(:,2) = R.*(cos(thetasim_interp(:,1)));

% Compare predictions to data for the score
RMSE=sqrt(mean((ysim-yexp(pureTriStarts:pureTriEnds,:)).^2,1));
f=sum(RMSE);

if nargout > 1
    % Compute the derivative of the objective
end

end

%%
function thetasim_interp = getCoordinates(params, data, fixed, st)
time  = data.time;
pureTriStarts = data.pureTriStarts(1);
pureTriEnds = data.pureTriEnds(1);

% Current parameters
R      = params(1);
AP      = params(2);
omega   = params(3);

% Calculate alpha
xzPosInit = data.com(pureTriStarts,:);
alpha = -asin((AP-xzPosInit(1))/R);

% Fixed quantities
g = fixed.g;
% Specify the ODE
ode = @(t,y) govEqu(y, g, R);
% Solve ODE over time interval
[tsim,thetasim] = ode45(ode,[time(pureTriStarts) time(pureTriEnds)],[alpha omega]); % first square bracket is the time range, second is initial condition
% Interpolate the simulated solution to predicted values at time points

thetasim_interp(:,1) = interp1(tsim,thetasim(:,1),time(pureTriStarts:pureTriEnds)); %For angle
thetasim_interp(:,2) = interp1(tsim,thetasim(:,2),time(pureTriStarts:pureTriEnds)); %For anglure velocity
end

%%
function dydt = govEqu(y, g, R)
%For SLIP
%y(1) = angle, y(2) = anglur speed
dydt =[y(2);
      g*sin(y(1))/R];
end