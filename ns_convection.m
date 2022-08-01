close all
clear
clc

%% Set up mesh and flow parameters
N = 41; % No of cols
M = 41; % No of rows
% N = 101;
% M = 101;
DX = 1/(N-1);
DY = 1/(M-1);
Pr = 0.7;
% Ra = 3*10^3;
Ra = 2*10^4;

%% Set up iteration parameters
% MaxItr = 1000000;
Itr = 1;
EPS1 = 10^(-3);
EPS2 = 10^(-3);
EPS3 = 10^(-3);
DT = 0.000001;

%% Set up initial flow
u = zeros(M,N);
v = zeros(M,N);

T = zeros(M,N);
T(1,:) = 1; % To meet boundary conditions (Only  2 to N-1?)
dTdX = zeros(M,N);
dTdY = zeros(M,N);
QT = zeros(M,N);

Psi = zeros(M,N);
SPsi = zeros(M,N);

Vor = zeros(M,N);
RVor = zeros(M,N);

%% Main  
while(1)

    RVor = Solve_RVor(N,M,Pr,Ra,Vor,RVor,T,u,v,DX,DY); % Update 2:N-1, 2:M-1
    Vor = Solve_Vor(N,M,Vor,RVor,DT); % Update 2:N-1, 2:M-1

    QT = Solve_QT(N,M,T,QT,u,v,DX,DY); % Update 2:N-1, 2:M-1
    T = Solve_Q(N,M,T,QT,DT); % Update 2:N-1, 2:M-1

    SPsi = Solve_SPsi(N,M,Vor,Psi,SPsi,DX,DY); % Update 3:N-2, 3:M-2
    Psi = Solve_Psi(N,M,DX,DY,SPsi,Psi); % Update 3:N-2, 3:M-2

    Psi = Update_BC_Psi(N,M,Psi); % Update 2nd and N-1/M-1th rows/columns (1st and N/Mth rows/columns are unchanged as they are bounday conditions)
    Vor = Update_BC_Vor(N,M,DX,DY,Vor,Psi); % Update 1st and N/Mth rows/columns
    T = Update_BC_T(N,M,T); % Update top-most and bottom-most rows (left and right most unchanges as boundary conditions)

    [u,v,dTdX,dTdY] = Compute_Vel_and_DT(N,M,DX,DY,Psi,T,u,v,dTdX,dTdY);

    [ErrorVor, ErrorT, ErrorPsi] = ComputeErrors(N,M,RVor,QT,SPsi);
    
    if mod(Itr,100) == 0
        disp(['Ietration: ', num2str(Itr), ': ErrorVor: ', num2str(ErrorVor), ', ErrorT: ', num2str(ErrorT), ', ErrorPsi: ', num2str(ErrorPsi)])
    end
    
    if (ErrorVor<=EPS1) && (ErrorT<=EPS2) && (ErrorPsi<=EPS3)
        disp("Convergence reached")
        break
    end
%     if Itr > MaxItr
%         disp("Hit max iteration bound")
%         break
%     end

    Itr = Itr + 1;
end

%% Simulation results
figure();
contourf(linspace(0,1,N),linspace(0,1,M),T);
colorbar;
xlabel('x')
ylabel('y')
title("Numerical solution T (grid: "+N+"by"+M+" , Ra: "+Ra+")");

figure();
contourf(linspace(0,1,N),linspace(0,1,M),u);
colorbar;
xlabel('x')
ylabel('y')
title("Numerical solution u (grid: "+N+"by"+M+" , Ra: "+Ra+")");

figure();
contourf(linspace(0,1,N),linspace(0,1,M),v);
colorbar;
xlabel('x')
ylabel('y')
title("Numerical solution v (grid: "+N+"by"+M+" , Ra: "+Ra+")");

figure();
contourf(linspace(0,1,N),linspace(0,1,M),Psi);
colorbar;
xlabel('x')
ylabel('y')
title("Numerical solution Psi (grid: "+N+"by"+M+" , Ra: "+Ra+")");

[umax, vmax, Nu0, Nu05, Numax, Numax_idx] = PerformTasks(N,M,u,v,dTdX);
disp([num2str(umax),';',num2str(vmax),';',num2str(Nu0),';',num2str(Nu05),';',num2str(Numax),' at j=',num2str(Numax_idx)])
%% Functions 
function RVor = Solve_RVor(N,M,Pr,Ra,Vor,RVor,T,u,v,DX,DY) 
    for j = 2:M-1
        for i = 2:N-1
            DVorX2 = (Vor(j,i+1)-2*Vor(j,i)+Vor(j,i-1))/(DX^2);
            DVorY2 = (Vor(j+1,i)-2*Vor(j,i)+Vor(j-1,i))/(DY^2);
            DVorX1 = u(j,i)*(Vor(j,i+1)-Vor(j,i-1))/(2*DX);
            DVorY1 = v(j,i)*(Vor(j+1,i)-Vor(j-1,i))/(2*DY);
            DTX1 = T(j,i+1)-T(j,i-1)/(2*DY);
            RVor(j,i) = Pr*(DVorX2+DVorY2)-DVorX1-DVorY1-Pr*Ra*DTX1;
        end
    end
end

function Vor = Solve_Vor(N,M,Vor,RVor,DT)
    for j = 2:M-1
        for i = 2:N-1
            Vor(j,i)=Vor(j,i)+DT*RVor(j,i);  
        end
    end
end

function QT = Solve_QT(N,M,T,QT,u,v,DX,DY)
    for j = 2:M-1
        for i = 2:N-1
            DTX2 = (T(j,i+1)-2*T(j,i)+T(j,i-1))/(DX^2);
            DTY2 = (T(j+1,i)-2*T(j,i)+T(j-1,i))/(DY^2);
            DTX1 = u(j,i)*(T(j,i+1)-T(j,i-1))/(2*DX);
            DTY1 = v(j,i)*(T(j+1,i)-T(j-1,i))/(2*DY);
            QT(j,i) = DTX2+DTY2-DTX1-DTY1;
        end
    end
end

function T = Solve_Q(N,M,T,QT,DT)
    for j = 2:M-1
        for i = 2:N-1
            T(j,i)=T(j,i)+DT*QT(j,i); 
        end
    end
end

function SPsi = Solve_SPsi(N,M,Vor,Psi,SPsi,DX,DY)
    for j = 3:M-2
        for i = 3:N-2
            SPsi(j,i)=Vor(j,i)-...
                      (Psi(j,i+1)-2*Psi(j,i)+Psi(j,i-1))/DX^2-...
                      (Psi(j+1,i)-2*Psi(j,i)+Psi(j-1,i))/DY^2;
        end
    end
end

function Psi = Solve_Psi(N,M,DX,DY,SPsi,Psi) 
    beta = 0.4; % Relaxation factor
    b_w = 1/(DX^2);
    b_s = 1/(DY^2);
    b_p = -2*(b_w+b_s);
    for j = 3:M-2
        for i = 3:N-2
            Psi(j,i)=Psi(j,i)+beta*SPsi(j,i)/b_p;
        end
    end
end

function Psi = Update_BC_Psi(N,M,Psi) % Update 2nd and N-1/M-1th rows/columns
    for j=2:M-1
        Psi(j,2) = 0.25*Psi(j,3);
        Psi(j,N-1) = 0.25*Psi(j,N-2);
    end
    
    for i=2:N-1
        Psi(2,i) = 0.25*Psi(3,i);
        Psi(M-1,i) = 0.25*Psi(M-2,i);
    end
end

function Vor = Update_BC_Vor(N,M,DX,DY,Vor,Psi) % Update 1st and N/Mth rows/columns
    for j=1:M
        Vor(j,1) = 3*Psi(j,2)/(DX^2) - 0.5*Vor(j,2);
        Vor(j,N) = 3*Psi(j,N-1)/(DX^2) - 0.5*Vor(j,N-1);
    end
    for i=1:N
        Vor(1,i) = 3*Psi(2,i)/(DY^2) - 0.5*Vor(2,i);
        Vor(N,i) = 3*Psi(M-1,i)/(DY^2) - 0.5*Vor(M-1,i);
    end
end

function T = Update_BC_T(N,M,T) % Update top-most and bottom-most rows
    for i=1:N
        T(i,1) = (4*T(i,2)-T(i,3))/3; 
        T(i,M) = (4*T(M-1,i)-T(M-2,i))/3;
    end
end

function [u,v,dTdX,dTdY] = Compute_Vel_and_DT(N,M,DX,DY,Psi,T,u,v,dTdX,dTdY)
    % Set physical boundary conditions for u, v
    % Set physical boundary conditions for dTdX, dTdY(top and bottom only)
    for j = 1:M
        u(j,1)=0;
        u(j,N)=0;
        v(j,1)=0;
        v(j,N)=0;
    end
    
    for i = 1:N
        u(1,i)=0;
        u(M,i)=0;
        v(1,i)=0;
        v(M,i)=0;
        dTdX(1,i)=0;
        dTdX(M,i)=0;
        dTdY(1,i)=0;
        dTdY(M,i)=0;
    end
    
    % Compute internal points
    for j = 2:M-1
        for i = 2:N-1
            u(j,i)=0.5*(Psi(j+1,i)-Psi(j-1,i))/DY;
            v(j,i)=0.5*(Psi(j,i-1)-Psi(j,i+1))/DX;
            dTdX(j,i)=0.5*(T(j,i+1)-T(j,i-1))/DX;
            dTdY(j,i)=0.5*(T(j+1,i)-T(j-1,i))/DY;
        end
    end
    
    % Compute dTdX at left and right most column
    for j = 1:M
        dTdX(j,1)=0.5*(4*T(j,2)-T(j,3))/DX;
        dTdX(j,N)=0.5*(-4*T(j,N-1)+T(j,N-2))/DX;
    end

end

function [ErrorVor, ErrorT, ErrorPsi] = ComputeErrors(N,M,RVor,QT,SPsi)
    ErrorVor = max(max(abs(RVor(2:M-1,2:N-1))));
    ErrorT = max(max(abs(QT(2:M-1,2:N-1))));
    ErrorPsi = max(max(abs(SPsi(3:M-2,3:N-2))));
end

function [umax, vmax, Nu0, Nu05, Numax, Numax_idx] = PerformTasks(N,M,u,v,dTdX)
    umax = max(u(:,round(N/2)));
    vmax = max(v(round(M/2),:));
    Nu0 = mean(dTdX(:,1));
    Nu05 = mean(dTdX(:,round(N/2)));
    [Numax,Numax_idx] = max(dTdX(:,1));    
end
