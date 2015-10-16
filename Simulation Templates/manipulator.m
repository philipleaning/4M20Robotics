%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%          3-Link Manipulator                                   %%%%%%
%%%%%%          M20 Robotics, coursework template                    %%%%%%
%%%%%%          University of Cambridge                              %%%%%%
%%%%%%          Michaelmas 2015                                      %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main function                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function manipulator
%% Blank slate
clear all
close all
clc

%% Parameters
l1 = 1;                 %length bottom link [m]
l2 = 1;                 %length middle link [m]
l3 = 0.5;               %length upper link [m]
lg = 0.2;               %length gripper [m]
wg = 0.2;               %width gripper [m]

%% Initial conditions
phi1 = 0;               %angle bottom link [rad]
phi2 = pi/2;            %angle middle link [rad]
phi3 = pi/2;            %angle upper link [rad]

phiVec = [phi1;phi2;phi3];
lVec = [l1,l2,l3];

%% Initialise animation
lineObj = animInit;
animation([phiVec,phiVec],lVec,lineObj,[-1,0.4]);

%% Run main routine
for i = 1:5

    pDes = [ginput(1)';0;1];    %get user defined target              
    phiVecNew = JacInv(phiVec,lVec,pDes);   %inverse kinematics
    animation([phiVec,phiVecNew],lVec,lineObj,pDes) %animation
    phiVec = phiVecNew;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Animation initialisation
function lineObj = animInit()
    
    lWin = l1+l2+l3+0.5;    %window size
    bOffs = 0.1;            %base offset
    
    figure(1)
    hold on
    % Ground plot
    plot([-lWin lWin],[-2*bOffs -2*bOffs],'k','LineWidth',3)
    plot(0,-bOffs,'^k','LineWidth',5)
    for j = 1:30
        
        k = -lWin-1+j*2*lWin/20;
        plot([k;k+0.2],[0-2*bOffs;-0.2-2*bOffs],'k','LineWidth',3)
    end

    %target line object
    lineObj.m1 = line(0,0,'color','r','LineWidth',5,'Marker','o');

    %link line objects
    lineObj.h1 = line(0,0,'color','k','LineWidth',2);
    lineObj.h2 = line(0,0,'color','k','LineWidth',2);
    lineObj.h3 = line(0,0,'color','k','LineWidth',2);

    %joint line objects
    lineObj.d1 = line(0,0,'color','k','LineWidth',5,'Marker','o');
    lineObj.d2 = line(0,0,'color','k','LineWidth',5,'Marker','o');

    %gripper line objects
    lineObj.g1 = line(0,0,'color','k','LineWidth',2);
    lineObj.g2 = line(0,0,'color','k','LineWidth',2);
    lineObj.g3 = line(0,0,'color','k','LineWidth',2);
    
    axis equal
    axis([-3 3 -0.4 3])
end

%% Animation
function animation(phiMat,lVec,lineObj,PDes)
    
    PhiVec = phiMat(:,2);
    PhiVecOld = phiMat(:,1);
    
    [A10,A21,A32] = HomCoord(PhiVec,lVec); %calculate transformation matrices
    
    r11 = A10*[0;l1;0;1];   %position link 1
    r22 = A10*A21*[0;l2;0;1];   %position link 2
    r33 = A10*A21*A32*[0;l3;0;1];   %position link 3

    rg11 = A10*A21*A32*[wg/2;l3;0;1];   %position gripper base right
    rg22 = A10*A21*A32*[wg/2;lg+l3;0;1];    %position gripper right
    rg11m = A10*A21*A32*[-wg/2;l3;0;1]; %position gripper base left
    rg22m = A10*A21*A32*[-wg/2;lg+l3;0;1];  %position gripper left
    
    %set target point
    set(lineObj.m1,'xdata',PDes(1),'ydata',PDes(2))
    
    %number of animated frames from start to end position
    n = 30;
    
    if PhiVecOld ~= PhiVec
        
        for k = 1:n

            %avoid multiple rotations
            PhiVec = mod(PhiVec,2*pi);  
            PhiVecOld = mod(PhiVecOld,2*pi);
            dPhi = PhiVec-PhiVecOld;
            dC = abs(dPhi)>pi;
            dPhi = dPhi - dPhi.*dC - sign(dPhi.*dC).*(2*pi*dC-abs(dPhi.*dC));

            %calculate intermediate positions
            phiVecT = PhiVecOld + dPhi*k/n;
            [A10T,A21T,A32T] = HomCoord(phiVecT,lVec);
            r1T = A10T*[0;l1;0;1];
            r2T = A10T*A21T*[0;l2;0;1];
            r3T = A10T*A21T*A32T*[0;l3;0;1];

            rg1T = A10T*A21T*A32T*[wg/2;l3;0;1];
            rg2T = A10T*A21T*A32T*[wg/2;lg+l3;0;1];
            rg1Tm = A10T*A21T*A32T*[-wg/2;l3;0;1];
            rg2Tm = A10T*A21T*A32T*[-wg/2;lg+l3;0;1];

            %set intermediate positions and pause
            set(lineObj.h1,'xdata',[0 r1T(1)],'ydata',[0 r1T(2)])
            set(lineObj.h2,'xdata',[r1T(1) r2T(1)],'ydata',[r1T(2) r2T(2)])
            set(lineObj.h3,'xdata',[r2T(1) r3T(1)],'ydata',[r2T(2) r3T(2)])

            set(lineObj.d1,'xdata',r1T(1),'ydata',r1T(2))
            set(lineObj.d2,'xdata',r2T(1),'ydata',r2T(2))

            set(lineObj.g1,'xdata',[rg1Tm(1) rg1T(1)],'ydata',[rg1Tm(2) rg1T(2)])
            set(lineObj.g2,'xdata',[rg1T(1) rg2T(1)],'ydata',[rg1T(2) rg2T(2)])
            set(lineObj.g3,'xdata',[rg1Tm(1) rg2Tm(1)],'ydata',[rg1Tm(2) rg2Tm(2)])
            pause(1/n)
        end
    end
    
    %set target position and pause
    set(lineObj.h1,'xdata',[0 r11(1)],'ydata',[0 r11(2)])
    set(lineObj.h2,'xdata',[r11(1) r22(1)],'ydata',[r11(2) r22(2)])
    set(lineObj.h3,'xdata',[r22(1) r33(1)],'ydata',[r22(2) r33(2)])

    set(lineObj.d1,'xdata',r11(1),'ydata',r11(2))
    set(lineObj.d2,'xdata',r22(1),'ydata',r22(2))

    set(lineObj.g1,'xdata',[rg11m(1) rg11(1)],'ydata',[rg11m(2) rg11(2)])
    set(lineObj.g2,'xdata',[rg11(1) rg22(1)],'ydata',[rg11(2) rg22(2)])
    set(lineObj.g3,'xdata',[rg11m(1) rg22m(1)],'ydata',[rg11m(2) rg22m(2)])

    drawnow
end

%% Homogeneous coordinate transform
function [A10,A21,A32] = HomCoord(phiVec,lVec)
    
    Phi1 = phiVec(1);
    Phi2 = phiVec(2);
    Phi3 = phiVec(3);

    L1 = lVec(1);
    L2 = lVec(2);

    A10 = [cos(Phi1) , -sin(Phi1) , 0 , 0; ...
           sin(Phi1) ,  cos(Phi1) , 0 , 0; ...
           0          , 0           , 1 , 0 ; ...
           0          , 0           , 0 , 1];

    A21 = [cos(Phi2) , -sin(Phi2) , 0 , 0; ...
           sin(Phi2) ,  cos(Phi2) , 0 , L1; ...
           0          , 0           , 1 , 0; ...
           0          , 0           , 0 , 1];

    A32 = [cos(Phi3) , -sin(Phi3) , 0 , 0; ...
           sin(Phi3) ,  cos(Phi3) , 0 , L2; ...
           0          , 0           , 1 , 0; ...
           0          , 0           , 0 , 1];

end

%% Jacobian inverse method (inverse kinematics)
function PhiVec = JacInv(PhiVec,LVec,PDes)
    
L3 = LVec(3)+0.1;   %offset target
j = 0;
[A10,A21,A32] = HomCoord(PhiVec,LVec);
rCurr = A10*A21*A32*[0;l3;0;1]; %current end effector position

pErr = norm(PDes-rCurr);    %current position error
errThresh = 0.001;  %error threshold

    while pErr > errThresh
        
        J = Jac(PhiVec,LVec);      
        dx =  pinv(J)*(PDes-rCurr); %pinv Moore-Penrose pseudoinverse
        PhiVec = PhiVec + dx;

        [A10,A21,A32] = HomCoord(PhiVec,LVec);

        rCurr = A10*A21*A32*[0;L3;0;1];
        pErr = norm(PDes-rCurr);

        j = j+1;

        if j>1000   %interrupt if position not below threshold after 1000 iterations

            %run optimisation routine to minimise position error
            lambda = fminsearch(@(x) JacErr(x,PhiVec,PDes,LVec),[0;0;0]);
            PhiVec = PhiVec + lambda;
            break
        end
    end
end

%% Calculate Jacobian
function J = Jac(phiVec,lVec)
    
phi1 = phiVec(1);
phi2 = phiVec(2);
phi3 = phiVec(3);

l1 = lVec(1);
l2 = lVec(2);
l3 = lVec(3);

% Jacobian
J = [ - l2*(cos(phi1)*cos(phi2) - sin(phi1)*sin(phi2)) - l1*cos(phi1) - ...
    l3*(cos(phi3)*(cos(phi1)*cos(phi2) - sin(phi1)*sin(phi2)) -...
    sin(phi3)*(cos(phi1)*sin(phi2) + cos(phi2)*sin(phi1))), -...
    l2*(cos(phi1)*cos(phi2) - sin(phi1)*sin(phi2)) - ...
    l3*(cos(phi3)*(cos(phi1)*cos(phi2) - sin(phi1)*sin(phi2)) -...
    sin(phi3)*(cos(phi1)*sin(phi2) + cos(phi2)*sin(phi1))), ...
    -l3*(cos(phi3)*(cos(phi1)*cos(phi2) - sin(phi1)*sin(phi2)) - ...
    sin(phi3)*(cos(phi1)*sin(phi2) + cos(phi2)*sin(phi1)));
    - l2*(cos(phi1)*sin(phi2) + cos(phi2)*sin(phi1)) - l1*sin(phi1) - ...
    l3*(cos(phi3)*(cos(phi1)*sin(phi2) + cos(phi2)*sin(phi1)) + ...
    sin(phi3)*(cos(phi1)*cos(phi2) - sin(phi1)*sin(phi2))), -...
    l2*(cos(phi1)*sin(phi2) + cos(phi2)*sin(phi1)) - ...
    l3*(cos(phi3)*(cos(phi1)*sin(phi2) + cos(phi2)*sin(phi1)) + ...
    sin(phi3)*(cos(phi1)*cos(phi2) - sin(phi1)*sin(phi2))), -...
    l3*(cos(phi3)*(cos(phi1)*sin(phi2) + cos(phi2)*sin(phi1)) +...
    sin(phi3)*(cos(phi1)*cos(phi2) - sin(phi1)*sin(phi2)));...
    0,    0,  0;...
    0,    0,  0];
end

%% Objective function error minimisation
function pErr = JacErr(lambda,phiVec,pDes,lVec)

l3 = lVec(3);
phi = phiVec + lambda;
[A10,A21,A32] = HomCoord(phi,lVec);
rCurr = A10*A21*A32*[0;l3;0;1];
pErr = norm(pDes-rCurr); 
end

end

