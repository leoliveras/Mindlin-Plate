%................................................................
% MATLAB codes for Finite Element Analysis
% problem19vibrations.m
% Mindlin plate in free vibrations
% antonio ferreira 2008

% clear memory
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
disp('Numero de GDL:')
disp(GDof)
disp('Numero de GDL prescritos:') % 'cccc' --> (3 GDL) *[4*(NElem+1)-4] Nos nas 4 bordas
disp(length(prescribedDof))
disp('Numero de nos ativos:')     % activeDof = prescribedDof - prescribedDof
disp(length(activeDof))


% load
F0 = -1;
freq = 2;
omega = 2*pi*freq;

% time
dt = 0.0001;
tmax = 2;
time = 0:dt:tmax;
nt = length(time);

% Newmark
gamma = 1/2;
betaN = 1/4;
ndof = length(activeDof);

U = zeros(ndof,nt);
Vt = zeros(ndof,nt);
A = zeros(ndof,nt);

K = stiffness(activeDof,activeDof);
M = mass(activeDof,activeDof);
C = dump(activeDof,activeDof);
Keff = K + gamma/(betaN*dt)*C + M/(betaN*dt^2);

for n = 1:nt-1

    % força física no instante atual
    P = F0*sin(omega*time(n));
    
    force = formForceVectorMindlinQ4( ...
        GDof, numberElements, ...
        elementNodes, numberNodes, ...
        nodeCoordinates, P);
    
    F = force(activeDof); % sum[F(:,1)] = F0 *sin

    Feff = F ...
        + M*( U(:,n)/(betaN*dt^2) ...
        + Vt(:,n)/(betaN*dt) ...
        + (1/(2*betaN)-1)*A(:,n) ) ...
        + C*( gamma/(betaN*dt)*U(:,n) ...
        + (gamma/betaN-1)*Vt(:,n) ...
        + dt*(gamma/(2*betaN)-1)*A(:,n) );

    U(:,n+1) = Keff\Feff;

    A(:,n+1) = (U(:,n+1)-U(:,n))/(betaN*dt^2) ...
             - Vt(:,n)/(betaN*dt) ...
             - (1/(2*betaN)-1)*A(:,n);

    Vt(:,n+1) = Vt(:,n) ...
              + dt*((1-gamma)*A(:,n)+gamma*A(:,n+1));

end


Uglobal = zeros(GDof,nt);
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