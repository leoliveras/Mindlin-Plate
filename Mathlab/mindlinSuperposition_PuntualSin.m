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
axis off
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
% sinusoidal
F0 = -1;
freq = 2;
omega = 2*pi*freq;

%position
centralNode = ceil(numberNodes/2);
forceDofGlobal = centralNode;
forceDof = find(activeDof==forceDofGlobal);
% disp('No carregado:')
% disp(centralNode)
% disp('Coordenadas:')
% disp(nodeCoordinates(centralNode,:))
% disp('GDL global da carga:')
% disp(forceDofGlobal)
% disp('Indice no vetor ativo:')
% disp(forceDof)
% disp('GDL ativo correspondente:')
% disp(activeDof(forceDof))

M = mass(activeDof,activeDof);
K = stiffness(activeDof,activeDof);
C = dump(activeDof,activeDof);

% V : mode shape
% D : frequency
%
numberOfModes = length(activeDof);
[V, D] = eigs(K, M, numberOfModes, 'smallestabs');
D=diag(sqrt(D)*L*sqrt(rho/G)); % Extrair os autovalores (frequências angulares ao quadrado)
[D,ii]=sort(D);  % Ordenar as frequências da menor para a maior
D = D(1:numberOfModes); % Selecionar apenas os primeiros N modos
ii = ii(1:numberOfModes);
V = V(:, ii); % Reordenar a matriz de autovetores (Modos de Vibrar)

% Normalização Robusta (Gram-Schmidt + Mass Norm)
% Isso força a ortogonalidade mesmo em modos repetidos (placas quadradas)
for i = 1:numberOfModes
    for j = 1:i-1
        % Projeção para garantir ortogonalidade (Gram-Schmidt)
        V(:,i) = V(:,i) - (V(:,j)' * M * V(:,i)) * V(:,j);
    end
    % Normalização para garantir o "1" na diagonal
    V(:,i) = V(:,i) / sqrt(V(:,i)' * M * V(:,i));
end

% Matrix projection
Mmod = V'*M*V;
Kmod = V'*K*V;
Cmod = V'*C*V;

%disp(Mmod(1:5,1:5))
%disp(Kmod(1:5,1:5))

% time
dt = 0.00001;
tmax = 2;
time = 0:dt:tmax;
nt = length(time);

Q = zeros(numberOfModes,nt);
Qd = zeros(numberOfModes,nt);
Qa = zeros(numberOfModes,nt);

% Newmark modal superposition
gamma = 1/2;
betaN = 1/4;

% força física no instante atual
Fphys = zeros(length(activeDof),1);

for n = 1:nt-1

    Fphys(forceDof,1) = F0*sin(omega*time(n));
    % projeção modal (calculada uma vez)
    Fmod = V'*Fphys;

    % integra cada modo
    for i = 1:numberOfModes

        m = Mmod(i,i);
        c = Cmod(i,i);
        k = Kmod(i,i);

        Keff = k + gamma/(betaN*dt)*c + m/(betaN*dt^2);

        f = Fmod(i);

        feff = f ...
            + m*( Q(i,n)/(betaN*dt^2) ...
            + Qd(i,n)/(betaN*dt) ...
            + (1/(2*betaN)-1)*Qa(i,n) ) ...
            + c*( gamma/(betaN*dt)*Q(i,n) ...
            + (gamma/betaN-1)*Qd(i,n) ...
            + dt*(gamma/(2*betaN)-1)*Qa(i,n) );

        Q(i,n+1) = feff/Keff;

        Qa(i,n+1) = (Q(i,n+1)-Q(i,n))/(betaN*dt^2) ...
                  - Qd(i,n)/(betaN*dt) ...
                  - (1/(2*betaN)-1)*Qa(i,n);

        Qd(i,n+1) = Qd(i,n) ...
                  + dt*((1-gamma)*Qa(i,n)+gamma*Qa(i,n+1));

    end
end

U = V*Q;

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