function [c,dcx, U] = Cal_c_ce_whole(nelx, nely, KE, xPhys, Emin, Emax, penal, freedofs, F)

nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(Emax-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK);
K = (K + K') / 2;
% fp = 2*([1 :  nelx+1] - 1)*(nely+1)+2;
% ff = zeros(1, length(fp))-0.05; ff(1) = -0.025;%ff(end)=-0.2;
% F = sparse( 2*(nelx*(nely+1)+nely/2+1), 1, -1, 2*(nely+1)*(nelx+1), 1);
% F = sparse( 2*((nelx+1)*(nely+1)), 1, -1, 2*(nely+1)*(nelx+1), 1);
%%
U = zeros(2*(nely+1)*(nelx+1), 1);
% fixeddofs = 1 : 2*(nely+1);
% tss = 2*(nelx*(nely+1)+1)-1 : 2 : 2 * (nelx+1)*(nely+1);
% fixeddofs = unique([fixeddofs, tss]);
% fixeddofs = 1 : 2*(nely+1);
% alldofs = 1 : 2*(nely+1)*(nelx+1); 
% freedofs = setdiff(alldofs, fixeddofs);
U(freedofs) = K(freedofs, freedofs)\F(freedofs);
ce = zeros(nely, nelx); 
ce(1 : nely, 1 : nelx) = reshape(sum((U(edofMat)*KE).*U(edofMat),2), nely, nelx);
c = sum(sum((Emin+reshape(xPhys(:), nely, nelx).^penal*(Emax-Emin)).*ce));
dcx = -penal*(Emax-Emin)*reshape(xPhys(:), nely, nelx).^(penal-1).*ce;

end