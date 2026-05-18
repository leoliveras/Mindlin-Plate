function [plotHandle] = getMeshNumeration(xx, yy)
    % Abre uma nova janela de figura
    figure;
    
    % Plota os nós e guarda o "handle" (identificador) do plot
    plotHandle = plot(xx, yy, 'ro', 'MarkerFaceColor', 'r');
    axis equal;
    grid on;
    hold on; % Garante que o texto seja desenhado sobre o plot

    % Opção A: Usando o tamanho do vetor xx para o loop
    numNodes = length(xx);
    for i = 1:numNodes
        % Adicionamos um pequeno deslocamento para o texto não ficar em cima do ponto
        text(xx(i), yy(i), [' ' num2str(i)], 'VerticalAlignment', 'bottom');
    end
end