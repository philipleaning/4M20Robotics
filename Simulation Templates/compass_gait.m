%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%          Compass gait walker                                  %%%%%%
%%%%%%          4M20 Robotics, coursework template                   %%%%%%
%%%%%%          University of Cambridge                              %%%%%%
%%%%%%          Michaelmas 2015                                      %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main function                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compass_gait
%% Blank slate
clear all
close all
clc

%% Parameters
g =9.8;
mh = 10;                %body mass [kg] 		suggested: 10
m = 5;                  %leg masses [kg]		suggested: 5
a = 0.5;                %'shin' length [m]      suggested: 0.5
b = 0.5;                %'thigh' length [m]     suggested: 0.5
gamma = 3;              %pitch of ramp [deg]	suggested: 3

gamma = gamma*pi/180;
l = a + b;

%% Initial Conditions
theta(1,1) = 0;     	%swing leg angle i/c [rad]              suggested: 0
theta(2,1) = 0;     	%stance leg angle i/c [rad]             suggested: 0
theta_dot(1,1) = 2;     %swing leg angular velocity i/c [rad/s]	suggested: 2
theta_dot(2,1) = -0.4;	%swing leg angular velocity i/c [rad/s]	suggested: -0.4

%% State variable
q = [theta; theta_dot];
qStore = [];

%% Simulation settings
number_steps = 20;  	%number of steps
zero_crossing = odeset('Events',@leg_switch);

%% Figure
figure(1)
title('Phase Portrait: Leg 1')
xlabel('theta [rad]')
ylabel('dtheta/dt [rad/s]')
axis([-0.4 0.5 -2.4 2.4])
hold on

%% Simulation
for i = 1:number_steps
	[Step(i).t,Step(i).q] = ode23(@EoM, [0,1], q', zero_crossing);
	
	alpha = abs((Step(i).q(end,1) - Step(i).q(end,2)))/2;
	
	q = Step(i).q(end,:)';
	
	W = W_matrix(alpha);
	
	q = W*q;
    
    qStore(end+1:end+length(Step(i).t),:) = [Step(i).t,Step(i).q,i*ones(length(Step(i).t),1)];

	if mod(i,2) == 0
		plot(Step(i).q(:,2),Step(i).q(:,4),'color',[0 0 0])
		line([Step(i).q(end,2) q(1)],[Step(i).q(end,4) q(3)],'color',[0 0 0])
	else
		plot(Step(i).q(:,1),Step(i).q(:,3),'color',[0 0 0])
		line([Step(i).q(end,1) q(2)],[Step(i).q(end,3) q(4)],'color',[0 0 0])
	end
	
end

Animation(qStore);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Equation of Motion
function output = EoM(~,q)
	
	theta(1,1) = q(1);
	theta(2,1) = q(2);
	theta_dot(1,1) = q(3);
	theta_dot(2,1) = q(4);
	
	H = H_matrix(theta);
	C = C_matrix(theta,theta_dot);
	G = G_matrix(theta);
	
	q_dot(1:2) = theta_dot;
	q_dot(3:4) =  -inv(H)*(C*theta_dot + G);
	
	output = q_dot';
end
%% H Matrix
function output = H_matrix(theta)
	output(1,1) = m*b^2;
	output(1,2) = -m*l*b*cos(theta(2)-theta(1));
	output(2,1) = -m*l*b*cos(theta(2)-theta(1));
	output(2,2) = (mh + m)*l^2 + m*a^2;
end
%% C Matrix
function output = C_matrix(theta,theta_dot)
	
	output(1,1) = 0;
	output(1,2) = m*l*b*sin(theta(2)-theta(1))*theta_dot(2);
	output(2,1) = -m*l*b*sin(theta(2)-theta(1))*theta_dot(1);
	output(2,2) = 0;
end
%% G Matrix
function output = G_matrix(theta)
	
	output(1,1) = m*b*g*sin(theta(1));
	output(2,1) = -(mh*l + m*a + m*l)*g*sin(theta(2));
end
%% Switch conditions
function [impact,terminate,direction] = leg_switch(~,q)

	if q(1) > 0
		impact = (q(1)+q(2)) + 2*gamma;
	else
		impact = 1;
	end
	terminate = 1;
	direction = -1;
end
%% W matrix
function output = W_matrix(alpha)
	
	J = [0 1; 1 0];
	
	Q_minus = Q_minus_matrix(alpha);
	Q_plus = Q_plus_matrix(alpha);
	
	H = (Q_plus)\Q_minus; % A\b == inv(A)*b
	
	output = [J zeros(2); zeros(2) H];
end
%% Q- matrix
function output = Q_minus_matrix(alpha)
	
	output(1,1) = -m*a*b;
	output(1,2) = -m*a*b + (mh*l^2 + 2*m*a*l)*cos(2*alpha);
	output(2,1) = 0;
	output(2,2) = -m*a*b;
end
%% Q+ matrix
function output = Q_plus_matrix(alpha)
	
	output(1,1) = m*b*(b - l*cos(2*alpha));
	output(1,2) = m*l*(l - b*cos(2*alpha)) + m*a^2 + mh*l^2;
	output(2,1) = m*b^2;
	output(2,2) = -m*b*l*cos(2*alpha);
end

function  Animation(qStore)
    %initialise animation
    figure(2)
  	l1 = line(0,0,'color','k','LineWidth',2);   %stance leg
    l2 = line(0,0,'color','k','LineWidth',2);   %swing leg
    gnd = line(0,0,'color','k','LineWidth',2);  %ground
    ml1 = line(0,0,'color','k','Marker','.','MarkerSize',m*10); %stance leg mass
    ml2 = line(0,0,'color','k','Marker','.','MarkerSize',m*10); %swing leg mass
    mlh = line(0,0,'color','k','Marker','.','MarkerSize',mh*10);%hip mass
    
    %sum up time and save in new column
    qStore(1,end+1) = qStore(1,1);
    for j = 2:length(qStore)
        if qStore(j,1) > qStore(j-1,1)
           qStore(j,end) = qStore(j-1,end) + abs(qStore(j,1)-qStore(j-1,1));
        else
           qStore(j,end) = qStore(j-1,end);
        end
    end
    
    %Initialise animation
    k = 2;
    rX = 0;
    rY = 0;
    RepSpeed = 1; %replay speed
    tic;
    while toc<qStore(end,end)/RepSpeed
        if toc > qStore(k,end)/RepSpeed
            %Store travelled distance after step ends
            if qStore(k,6) > qStore(k-1,6)     
                 rX = r_12(1);
                 rY = r_12(2);
            end
            r_f1 = [-(a+b)*sin(qStore(k,3));(a+b)*cos(qStore(k,3))];    %stance foot to hip vector
            r_O1 = [rX;rY] + r_f1;  %origin to hip vector
            r_f2 = [(a+b)*sin(qStore(k,2));-(a+b)*cos(qStore(k,2))];    %hip to swing foot vector
            r_12 = r_O1 + r_f2;  %origin to swing foot vector
            
            %set line objects
            set(l1,'xdata',[rX r_O1(1)],'ydata',[rY r_O1(2)])
            set(l2,'xdata',[r_O1(1) r_12(1)],'ydata',[r_O1(2) r_12(2)])
            set(ml1,'xdata',rX+r_f1(1)/norm(r_f1)*a,'ydata',rY + r_f1(2)/norm(r_f1)*a)
            set(ml2,'xdata',r_O1(1) + r_f2(1)/norm(r_f2)*b,'ydata',r_O1(2) + r_f2(2)/norm(r_f2)*b)
            set(mlh,'xdata',r_O1(1),'ydata',r_O1(2))
            gndX = (r_12(1)+rX)/2;
            gndY = -gndX*sin(gamma);
            set(gnd,'xdata',[gndX + cos(gamma), gndX-cos(gamma)],'ydata',[gndY - sin(gamma), gndY + sin(gamma)])
                       
            axis([gndX-1 gndX+1 -0.2+gndY 1.2+gndY])
            drawnow
            k = k+1;
        end
    end
end

end
