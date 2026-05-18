% clear memory
format long
clear; clc;

% materials
E = 10920; poisson = 0.30;
thickness=0.1;
I=thickness^3/12;
rho=1;
%
kapa=0.8601; % cccc / cccf case
%kapa=0.822; % scsc case
%kapa=5/6; % ssss case

% matrix B of kinematic compatibility
% bending part
B_bending=I*E/(1-poisson^2)*[1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];
% shear part
B_shear=kapa*thickness*E/2/(1+poisson)*eye(2);


% Mesh generation
L = 1;
numberElementsX=10;
numberElementsY=10;
numberElements=numberElementsX*numberElementsY;
[nodeCoordinates, elementNodes] = rectangularMesh(L,L,numberElementsX,numberElementsY);
xx=nodeCoordinates(:,1);
yy=nodeCoordinates(:,2);

drawingMesh(nodeCoordinates,elementNodes,'Q4','k-');
[ploteNumeration] = getMeshNumeration(xx,yy);
numberNodes=size(xx,1);

% GDof: global number of degrees of freedom
GDof=3*numberNodes;

% computation of the system stiffness, dumpi and mass matrices
[stiffness] = formStiffnessMatrixMindlinQ4(GDof,numberElements,...
                elementNodes,numberNodes,nodeCoordinates,B_shear,...
                B_bending,thickness,I);

[mass] = formMassMatrixMindlinQ4(GDof,numberElements,...
            elementNodes,numberNodes,nodeCoordinates,thickness,rho,I);

%dump = zeros(size(stiffness));
alpha = 0.1;
beta = 0.01;
[dump] = alpha*mass + beta*stiffness;

% % boundary conditions
[prescribedDof,activeDof,fixedNodeW]=EssentialBC('cccc',GDof,xx,yy,nodeCoordinates,numberNodes);
G=E/2.6;

% load
centralNode = ceil(numberNodes/2);
forceDofGlobal = centralNode;
forceDof = find(activeDof==forceDofGlobal);
F0 = -1;
freq = 2;
omega = 2*pi*freq;

% time convergence
dt = 1;
fn_max = 678.1079; %Hz
wn_max = 2*pi*fn_max;
dt_max= 2/wn_max;
if dt > dt_max
    %dt = dt_max;
    dt = 10^floor(log10(dt_max));
end

% time step
tmax = 2;
time = 0:dt:tmax;
nt = length(time);

% Central diferences
K = stiffness(activeDof,activeDof);
M = mass(activeDof,activeDof);
C = dump(activeDof,activeDof);

ndof = length(activeDof);
U = zeros(ndof,nt);

Uglobal = zeros(GDof,nt);
Vt = zeros(ndof,nt);
A = zeros(ndof,nt);

for n = 1:nt-1

    % load
    F = zeros(ndof,1);
    F(forceDof) = F0*sin(omega*time(n));
    
    A(:,1) = M\(F - C*Vt(:,1) - K*U(:,1));
    Uminus1 = U(:,1) - dt*Vt(:,1) + 0.5*dt^2*A(:,1);

    % Aux matrix
    A1 = M/dt^2 + C/(2*dt);
    A2 = K - 2*M/dt^2;
    A3 = M/dt^2 - C/(2*dt);

    if n == 1
        U(:,n+1) = A1 \ (F - A2*U(:,n) - A3*Uminus1);
    else
        U(:,n+1) = A1 \ (F - A2*U(:,n) - A3*U(:,n-1));
    end

    % velocidades e acelerações (opcional)
    if n > 1
        Vt(:,n) = (U(:,n+1)-U(:,n-1))/(2*dt);

        A(:,n) = (U(:,n+1)-2*U(:,n)+U(:,n-1))/dt^2;
    end
end

for n = 1:nt
    Uglobal(activeDof,n) = U(:,n);
end
wDof = 1:numberNodes;   %EssentialBC usa ordem global: primeiro todos os w, depois todos os θx, depois todos os θy
                        % ou seja: [ w1 … wn ∣  θx ∣  θy ]
disp('Deslocamento maximo:')
wmax = max(abs(Uglobal(wDof,:)),[],'all')
w = Uglobal(wDof,n);
figure
scatter(xx,yy,80,w,'filled')
colorbar
axis equal

X = reshape(xx,numberElementsY+1,numberElementsX+1);
Y = reshape(yy,numberElementsY+1,numberElementsX+1);


figure
for n = 1:500:nt
    clf
    w = Uglobal(wDof,n);
    W = reshape(w,numberElementsY+1,numberElementsX+1);
    surf(X,Y,W)
    xlabel('x')
    ylabel('y')
    zlabel('w')
    title(['t = ', num2str(time(n))])
    zlim([-wmax wmax])
    drawnow
end
 
% plot(time,U(forceDof,:))
% xlabel('Tempo')
% ylabel('w')
% grid on
99;