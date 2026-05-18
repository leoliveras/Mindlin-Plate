%................................................................
function displacements = solution(GDof, prescribed, stiffness,force)
    % function to find solution in terms of global displacements
    activeDof=setdiff([1:GDof]', [prescribed]);
    U=stiffness(activeDof,activeDof)\force(activeDof);
    displacements=zeros(GDof,1);
    displacements(activeDof)=U;
end