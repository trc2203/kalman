stream=fopen('sensor_data.txt','r'); % Load the simulated sensor data into the vectors time, mb, and sb
header=fscanf(stream,'%s',8);
i=1;
time=zeros(601,1);mb=zeros(601,3);sb=zeros(601,3); % Preallocate arrays for speed
while(feof(stream)==0)
    y=fscanf(stream,'%e,%e,%e,%e,%e,%e,%e\\n');
    time(i)=y(1);
    mb(i,1:3)=y(2:4);
    sb(i,1:3)=y(5:7);
    i=i+1;
end
fclose(stream); % Close the sensor data input file
stream2=fopen('kalman_filter_output.txt','w'); % Open the kalman filter output file and write the header information
fprintf(stream2,'time [s],q1,q2,q3,q4,wx,wy,wz\n');
N=length(time);
xhk=[1;0;0;0;0;0;0]; % Initial state vector input xhatk
Pk=eye(7);
R=deg2rad([eye(3),zeros(3);
    zeros(3),0.1*eye(3)]); % Error in measurements
dt=[0,0.1];
deltat=dt(2)-dt(1);
Q=eye(7); % Process noise 7x7 identity matrix
Qscale=1; % Process noise scaling factor
X=zeros(601,7);RPY=zeros(601,3);sumsquares=zeros(601,1); % Preallocate arrays for speed
for i=1:N
    yi=[mb(i,:),sb(i,:)]';
    [t,xdot]=ode45(@prop,dt,xhk); % Propagate the state
    xhkp1m=xdot(end,:)'; % Only want last line of xdot
    q1=xhkp1m(1);q2=xhkp1m(2);q3=xhkp1m(3);q4=xhkp1m(4);
    wx=xhkp1m(5);wy=xhkp1m(6);wz=xhkp1m(7);
    A=0.5*[0,wz,-wy,wx,q4,-q3,q2;
        -wz,0,wx,wy,q3,q4,-q1;
        wy,-wx,0,wz,-q2,q1,q4;
        -wx,-wy,-wz,0,-q1,-q2,-q3;
        zeros(3,7)]; % Linearise around the last estimate
    C=2*[-q3,q4,-q1,q2,0,0,0;
        -q4,-q3,-q2,-q1,0,0,0;
        q1,q2,-q3,-q4,0,0,0;
        -q1,q2,q3,-q4,0,0,0;
        -q2,-q1,q4,q3,0,0,0;
        -q3,-q4,-q1,-q2,0,0,0];
    Phi=eye(7)+A*deltat+deltat^2*A^2/factorial(2)+deltat^3*A^3/factorial(3);
    Pkp1m=Phi*Pk*Phi'+Q*Qscale; % Propagate the covariance (Qscale allows scaling process noise)
    K=Pkp1m*C'/(C*Pkp1m*C'+R); % Calculate the Kalman gain, / equivalent to inv() but more efficient
    h=[2*q2*q4-2*q1*q3;
        -2*q1*q4-2*q2*q3;
        q1^2+q2^2-q3^2-q4^2;
        -q1^2+q2^2+q3^2-q4^2;
        2*q3*q4-2*q1*q2;
        -2*q1*q3-2*q2*q4];
    xhkp1=xhkp1m+K*(yi-h); % Update the state estimate
    fprintf(stream2,'%g,%g,%g,%g,%g,%g,%g,%g\n',time(i),xhkp1'); % Output the results to file
    xhk=xhkp1; % Redefine x for next iteration
    Pk=Pkp1m-K*C*Pkp1m; % Update error covariance estimate
    X(i,:)=xhk; % Store xhat in an array to plot
    roll=atan2(2*(xhk(4)*xhk(1)+xhk(2)*xhk(3)),1-2*(xhk(1)^2+xhk(2)^2)); % Convert quaternions to Euler angles - Phi
    pitch=asin(2*(xhk(4)*xhk(2)-xhk(3)*xhk(1))); % Theta
    yaw=atan2(2*(xhk(4)*xhk(3)+xhk(1)*xhk(2)),1-2*(xhk(2)^2+xhk(3)^2)); % Psi
    RPY(i,:)=[roll,pitch,yaw];
    sumsquares(i)=X(i,1)^2+X(i,2)^2+X(i,3)^2+X(i,4)^2;
end
fclose(stream2); % Close the kalman filter output file
figure(1)
plot(time,X(:,1:4)) % Plot quaternions q1, q2, q3, q4
xlabel('Time (s)')
ylabel('Quaternions')
legend('q_1','q_2','q_3','q_4')
figure(2)
plot(time,X(:,5:7)) % Plot angular rates of spacecraft wx, wy, wz
xlabel('Time (s)')
ylabel('Estimated angular rates (rad s^{-1})')
legend('\omega_x','\omega_y','\omega_z','interpreter','latex')
figure(3)
plot(time,RPY) % Plot Euler angles phi, theta, psi
xlabel('Time (s)')
ylabel('Euler angles (rad)')
legend('\phi','\theta','\psi','interpreter','latex')
figure(4)
plot(time,sumsquares)
xlabel('Time (s)')
ylabel('Sum of the squares of the quaternions')
stream3=fopen('final.txt','w');
fprintf(stream3,'The final orientation of the satellite in quaternion [q1, q2, q3, q4] representation is q1=%g, q2=%g, q3=%g, q4=%g\n',X(end,1:4));
fprintf(stream3,'The final orientation of the satellite in Euler angle [roll, pitch, yaw] representation is phi=%g rad, theta=%g rad, psi=%g rad\n',RPY(end,:));
fprintf(stream3,'The final estimated angular rates of the satellite are w_x=%g rad s^-1, w_y=%g rad s^-1, w_z=%g rad s^-1\n',X(end,5:7));
fclose(stream3);