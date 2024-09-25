function [iKT,jKT,KET,TT,edofMatFT]=PreKTQ4(nelx, nely, coordinate)

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
    
%    KET2=KET2+W(i)*det(J0)*(N_xi*N_xi');
end

%% global stiffness
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVecT = reshape(nodenrs(1:end-1,1:end-1),nelx*nely,1);
edofMatFT = repmat(edofVecT,1,4) + ...
           repmat([1,nely+2,nely+1,0],nelx*nely,1);
% TF
iTT = reshape(edofMatFT,4*nelx*nely,1);
jTT = reshape(repmat([1:nelx*nely],4,1)',4*nelx*nely,1);
sTT = repmat(1/4,4*nelx*nely,1);
TT = sparse(iTT,jTT,sTT);

% KF
iKT = reshape(kron(edofMatFT,ones(4,1))',16*nelx*nely,1);
jKT = reshape(kron(edofMatFT,ones(1,4))',16*nelx*nely,1);

end