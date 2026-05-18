function drawEigenmodes2D(x,y,VVV,NN,N,D)
% Draw vibration modes for Mindlin plate
%
% x,y  : coordinate vectors
% VVV  : modal matrix (each column = one mode)
% NN   : total number of nodes
% N    : nodes per direction
% D    : natural frequencies

numberModes = size(VVV,2);

[X,Y] = meshgrid(x,y);

for mode = 1:numberModes

    figure

    % reshape modal displacement into 2D field
    W = reshape(VVV(:,mode),N,N);

    surf(X,Y,W)

    xlabel('x')
    ylabel('y')
    zlabel('w')

    title(sprintf(...
        'Mode %d   Frequency = %.4f', ...
        mode, D(mode)))

    shading interp
    axis equal
    view(45,30)

end

end