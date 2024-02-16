%	FLIGHT --  Generic 6-DOF Trim, Linear Model, and Flight Path Simulation
	clear
    clear all
	global GEAR CONHIS SPOIL u x V tuHis deluHis uInc TrimHist SMI MODEL RUNNING

    disp('** 6-DOF FLIGHT Simulation **')
    date
    
    MODEL   =   1;      % Aerodynamic model selection
                        % 0: Incompressible flow, high angle of attack
                        % 1: Compressible flow, low angle of attack    
	TRIM    = 	1;		% Trim flag (= 1 to calculate trim @ I.C.)
	LINEAR  = 	0;		% Linear model flag (= 1 to calculate and store F and G)
	SIMUL   =	1;		% Flight path flag (= 1 for nonlinear simulation)
	GEAR    = 	0;		% Landing gear DOWN (= 1) or UP (= 0)
	SPOIL   =	0;		% Symmetric Spoiler DEPLOYED (= 1) or CLOSED (= 0)
	CONHIS  =	1;		% Control history ON (= 1) or OFF (= 0)
	dF      = 	0;		% Flap setting, deg
	RUNNING =   0;      % internal flag, -
    
%	Initial Altitude (ft), Indicated Airspeed (kt), 
%   Dynamic Pressure (N/m^2), and True Airspeed (m/s

	hft         =   10000   % Altitude above Sea Level, ft
    VKIAS       =   150     % Indicated Airspeed, kt
    
        hm          =   hft * 0.3048    % Altitude above Sea Level, m
        VmsIAS      =   VKIAS * 0.5154  % Indicated Airspeed, m/s
    
        [airDens,airPres,temp,soundSpeed] = Atmos(hm)
    
        qBarSL  =   0.5*1.225*VmsIAS^2  % Dynamic Pressure at sea level, N/m^2
        V   =   sqrt(2*qBarSL/airDens);	% True Airspeed, TAS, m/s
        TASms   =   V

%	Alphabetical List of Initial Conditions

	alpha   =	0;      % Angle of attack, deg (relative to air mass)
	beta    =	0;      % Sideslip angle, deg (relative to air mass)
	dA      =	0;      % Aileron angle, deg
	dAS     =	0;      % Asymmetric spoiler angle, deg
	dE      =	0;      % Elevator angle, deg
	dR      =	0;      % Rudder angle, deg
	dS      = 	0;      % Stabilator setting, deg
	dT      = 	0;      % Throttle setting, % / 100
	hdot    =	0;      % Altitude rate, m/s
	p       =	0;      % Body-axis roll rate, deg/s
	phi     =	0;      % Body roll angle wrt earth, deg
	psi     =	0;      % Body yaw angle wrt earth, deg
	q       =	0;      % Body-axis pitch rate, deg/sec
	r       =	0;      % Body-axis yaw rate, deg/s
	SMI     =	0;      % Static margin increment due to center-of-mass variation from reference, %/100
	tf      =	100;    % Final time for simulation, sec
	ti      = 	0;      % Initial time for simulation, sec
	theta   =	alpha;  % Body pitch angle wrt earth, deg [theta = alpha if hdot = 0]
	xe      =	0;      % Initial longitudinal position, m
	ye      = 	0;      % Initial lateral position, m
	ze      = 	-hm;    % Initial vertical position, m [h: + up, z: + down]
    
    if MODEL == 0
        disp('Low-Mach, High-Alpha Model')
    else
        disp('High-Mach, Low Alpha Aerodynamic Model')
    end
    

%	Initial Conditions Depending on Prior Initial Conditions

	gamma	=	57.2957795 * atan(hdot / sqrt(V^2 - hdot^2))
						% Inertial Vertical Flight Path Angle, deg
	qbar	= 	0.5 * airDens * V^2	
						% Dynamic Pressure, N/m^2
	IAS		=	sqrt(2 * qbar / 1.225)
						% Indicated Air Speed, m/s
	Mach	= 	V / soundSpeed	
						% Mach Number
											
%	Initial State Perturbation (Test Inputs: m, m/s, rad, or rad/s)
	delx	=	[0;0;0
				0;0;0
                0;0;0
				0;0;0];
            
%	Initial Control Perturbation (Test Inputs: rad or 100%)			
	delu	=	[0;0;0;-0.234;0;0;0];
	
%	Control Perturbation History (Test Inputs: rad or 100%)
%   =======================================================
%   Each control effector represented by a column
%	Each row contains control increment delta-u(t) at time t:

	tuHis	=	[0 33 67 100];
	deluHis	=	[0 0 0 0 0 0 0
				0 0 0 0 0 0 0
			    0 0 0 0 0 0 0
			    0 0 0 0 0 0 0];
	uInc	=	[];

%	State Vector and Control Initialization, rad
	phir	=	phi * .01745329;
	thetar	=	theta * .01745329;
	psir	=	psi * .01745329;

	windb	=	WindField(-ze,phir,thetar,psir);
	alphar	=	alpha * .01745329;
	betar	=	beta * .01745329;

	x	=	[V * cos(alphar) * cos(betar) - windb(1)
			V * sin(betar) - windb(2)
			V * sin(alphar) * cos(betar) - windb(3)
			xe
			ye
			ze
			p * 0.01745329
			q * 0.01745329
			r * 0.01745329
			phir
			thetar
			psir];
	
	u	=	[dE * 0.01745329
			dA * 0.01745329
			dR * 0.01745329
			dT
			dAS * 0.01745329
			dF * 0.01745329
			dS * 0.01745329];

%	Trim Calculation (for Steady Level Flight at Initial V and h)

	if TRIM >= 1
		'Calculate TRIM Stabilator, Thrust, and Pitch Angle'
        OptParam        =   [];
        TrimHist        =   [];
		InitParam		=	[0.0369;0.1892;0.0986];
        
		[OptParam,J,ExitFlag,Output]  =	fminsearch('TrimCost',InitParam)
        TrimHist;
        Index=  [1:length(TrimHist)];
        TrimStabDeg     =   57.2957795*OptParam(1)
		TrimThrusPer    =   100*OptParam(2)
        TrimPitchDeg    =   57.2957795*OptParam(3)
        TrimAlphaDeg    =   TrimPitchDeg - gamma
        
%		Insert trim values in nominal control and state vectors
        'Incorporate trim values in control and state vectors'
		u	=	[u(1)
				u(2)
				u(3)
				OptParam(2)
				u(5)
				u(6)
				OptParam(1)]
		format long			
		x	=	[V * cos(OptParam(3))
				x(2)
				V * sin(OptParam(3))
				x(4)
				x(5)
				x(6)
				x(7)
				x(8)
				x(9)
				x(10)
				OptParam(3)
				x(12)]
		format short
    end
 

%	Flight Path Calculation

    if SIMUL >= 1
        RUNNING =   1;  
		tspan	=	[ti tf];
		xo		=	x + delx
		u		=	u + delu
        
        options =   odeset('Events',@event,'RelTol',1e-7,'AbsTol',1e-7);
		[t,x]	=	ode15s(@EoM,tspan,xo,options);
        
		kHis	=	length(t)

                VAirRel         =   [];
        vEarth          =   [];
        AlphaAR         =   [];
        BetaAR          =   [];
        windBody        =   [];
        airDensHis      =   [];
        soundSpeedHis   =   [];
        qbarHis         =   [];
        GammaHis        =   [];
        XiHis           =   [];
        
        for i = 1:kHis
            windb           =	WindField(-x(i,6),x(i,10),x(i,11),x(i,12));
            windBody        =   [windBody windb];
            [airDens,airPres,temp,soundSpeed] = Atmos(-x(i,6));
            airDensHis      =   [airDensHis airDens];
            soundSpeedHis   =   [soundSpeedHis soundSpeed];
        end
        
        vBody           =   [x(:,1) x(:,2) x(:,3)]';
        vBodyAir        =   vBody + windBody;

        for i = 1:kHis
            vE          =   DCM(x(i,10),x(i,11),x(i,12))' * [vBody(1,i);vBody(2,i);vBody(3,i)];
            VER         =   sqrt(vE(1)^2 + vE(2)^2 + vE(3)^2);
            VAR         =   sqrt(vBodyAir(1,i)^2 + vBodyAir(2,i)^2 + vBodyAir(3,i)^2);
            VARB        =   sqrt(vBodyAir(1,i)^2 + vBodyAir(3,i)^2);

            if vBodyAir(1,i) >= 0
                Alphar      =	asin(vBodyAir(3,i) / VARB);
            else
                Alphar      =	pi - asin(vBodyAir(3,i) / VARB);
            end
            
            AlphaAR     =   [AlphaAR Alphar];
            Betar       = 	asin(vBodyAir(2,i) / VAR);
            BetaAR      =   [BetaAR Betar];
            vEarth      =   [vEarth vE];
            Gammar      =   asin(-vEarth(3,i) / VER);
            GammaHis    =   [GammaHis Gammar];
            Xir         =   asin(vEarth(2,i) / sqrt((vEarth(1,i))^2 + (vEarth(2,i))^2));
            if vEarth(1,i) <= 0 && vEarth(2,i) <= 0
                Xir = -pi - Xir;
            end
            if vEarth(1,i) <= 0 && vEarth(2,i) >= 0
                Xir = pi - Xir;
            end
            
            XiHis       =   [XiHis Xir];
            VAirRel     =   [VAirRel VAR];
        end

        MachHis         =   VAirRel ./ soundSpeedHis;
        AlphaDegHis 	=	57.2957795 * AlphaAR;
        BetaDegHis      =	57.2957795 * BetaAR;
        qbarHis         =	0.5 * airDensHis .* VAirRel.*VAirRel;
        GammaDegHis     =   57.2957795 * GammaHis;
        XiDegHis        =   57.2957795 * XiHis;



% Data points
% dp= qbarHis'
numbers=MachHis'


% Time intervals
time_intervals =t;

% Plotting measured data
figure;
plot(time_intervals,numbers, 'b-', 'LineWidth', 1.5);
xlabel('Time(s)');
ylabel('velocity (Mach)');
grid on;

% Adding noise to the measured values
noise_stddev =0.003; % Standard deviation of the noise
noisy_numbers =  numbers+ noise_stddev * randn(size(numbers));

% Applying Kalman filter
% Assuming initial state and measurement noise covariance
initial_state = noisy_numbers(1);
measurement_noise_covariance = noise_stddev^2;

% Kalman filter initialization
x = initial_state; % Initial state estimate
P = measurement_noise_covariance; % Initial error covariance

filtered_values = zeros(size(noisy_numbers));
filtered_values(1) = x;

for i = 2:numel(noisy_numbers)
    % Prediction step
    x = x;
    P = P + measurement_noise_covariance;
    
    % Update step
    K = P / (P + measurement_noise_covariance);
    x = x + K * (noisy_numbers(i) - x);
    P = (1 - K) * P;
    
    filtered_values(i) = x;
end

% Time interval (assumed)
dt = 1; % You can adjust this value according to your assumption

% Generate time values
time =t;

% Perform least squares fitting
degree =6; % Quadratic curve

% Construct the Vandermonde matrix
A = ones(length(numbers), degree + 1);
for i = 1:degree
    A(:, i+1) = time.^i;
end

% Use the least squares method to find the coefficients
coefficients = A \ numbers;

% Generate the fitted curve
curve = A * coefficients;

% Plotting filtered data
hold on;
plot(time_intervals, filtered_values, 'r-', 'LineWidth', 1.5);

hold on

plot(t, curve, 'LineWidth', 2);
legend('Measured Data with Noise', 'Filtered Data','least square Data');
grid on
hold off

    end
