function [KET]=PreKTQ4_new(coordinate)

%% Element stiffness
quadorder=2;
[W,Q] = GAUSS(quadorder,2);

KET=zeros(4,4);
% KET2=zeros(4,4);
for i=1:size(Q,1)
    [N_xi, dNdxi] = NShape('Q4',Q(i,:));  % element shape functions
    
    J0 = dNdxi'*coordinate;           % element Jacobian matrix
    invJ0 = inv(J0);
    dNdx  = invJ0*dNdxi';             % derivatives of N w.r.t XY
    Bx=dNdx(1,:)'*dNdx(1,:);
    By=dNdx(2,:)'*dNdx(2,:);
    
    KET=KET+W(i)*det(J0)*(Bx+By);
    
   % KET2=KET2+W(i)*det(J0)*(N_xi*N_xi');
end

end