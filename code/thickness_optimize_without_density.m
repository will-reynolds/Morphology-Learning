clear
close all
tic
rmin = 6;
nelx = 120;
nely = 60;
nele = nely*nelx;
xPhys = ones(nely, nelx);
voidSet = [];
nStage = 20;
nloop = 1000;
volfrac = 0.6;
tolvol = sum(xPhys(:));
colormap(gray); imagesc(1-xPhys); axis equal; axis tight; axis off; drawnow;

%% CONNECTIVITY MATRIX / LAPLACE MATRIX
M = zeros(4, nely*nelx);
for x = 1 : nelx
    for y = 1 : nely
        i = nely*(x - 1) + y;
        M(:, i) = [(nely+1)*x + y; (nely+1)*(x-1) + y; (nely+1)*(x-1)+y + 1; (nely+1)*x + y + 1];
    end
end
%%
X = repmat(0.5:nelx+.5, nely+1, 1);
X = flipud(X);
Y = repmat([nely+.5:-1:0.5]', 1, nelx+1);
% Y = flipud(Y);
V = [X(:), Y(:)];
F = M';
VE = (V(F(:, 1), :) + V(F(:, 2), :) + V(F(:, 3), :) + V(F(:, 4), :))/4;

%% MATERIAL PROPERTIES
Emax = 1;
Emin = 1e-9; 
nu = 0.3;
E0 = zeros(1, nely*nelx) + Emin;
E0(voidSet) = 0;

%% INITIALIZE SVERAL PARAMETERS
penal = 3;      % stiffness penalty

%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

%% PREPARE FILTER DENSITY
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                if ismember(e1, voidSet) && e1 ~= e2
                    continue;
                elseif ismember(e2, voidSet) && ~ismember(e1, voidSet)
                    continue;
                end
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

%% CONNECTIVITY MATRIX / LAPLACE MATRIX
atria = nn_prepare(VE);
kk = 10;
[ind, dis] = nn_search(VE, atria, VE, kk);
II = zeros(nele*kk, 1);
JJ = zeros(nele*kk, 1);
SS = zeros(nele*kk, 1);
k = 0;
for i = 1 : size(VE, 1)
    [xx, yy] = find(dis(i, :) <= sqrt(1)+1.0e-6);
    t = setdiff(ind(i, yy), i);
    II(k+1:k+length(t)) = i;
    JJ(k+1:k+length(t)) = t;
    SS(k+1:k+length(t)) = 1;
    k = k + length(t);
end

[xx, yy] = find(II == 0);
II(xx) = []; JJ(xx) = []; SS(xx) = [];
L = sparse(II, JJ, SS);
a = max(sum(L, 2), 1.0e-19);
for i = 1 : size(L, 1)
    [xx, yy] = find(L(i, :) ~= 0);
    L(i, yy) = L(i, yy) / a(i);
end
L = sparse(L);
E = speye(size(L, 1));
LM = E - L;
LM = sparse(LM);

%% force
ff = sparse(2*(nely+1)*(nelx+1), 1);
ff(2*((nely+1)*nelx+1)) = -1;

%% fixed freedom of degree
fixedNode = 1:nely+1;
hold on
draw_points(V(fixedNode, :), 10, [1 1 0]);
cameratoolbar
[xx, yy] = find(xPhys == 1);
aa = [(yy-1)*(nely+1)+xx, (yy-1)*(nely+1)+xx+1, yy*(nely+1)+xx, yy*(nely+1)+xx+1];
aa = unique(aa(:));
bb = 1 : (nely+1)*(nelx+1);
aa = setdiff(bb, aa);
fixedNode = [fixedNode, aa];
fixeddofs = unique([2*fixedNode, 2*fixedNode-1]);
alldofs = 1 : 2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs, fixeddofs);

%% for heat equation
%% global stiffness
nodenrsh = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVecT = reshape(nodenrsh(1:end-1,1:end-1),nelx*nely,1);
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

%% INITIALIZE ITERATION
%%
low = 0;
upp = 0;
rou = 1;
loop = 0;
beta = 1;       % beta continuation
eta = 0.5;      % projection threshold, fixed at 0.5
% x = repmat(volfrac,nely,nelx);
x = rand(nely, nelx);
xTilde(:) = (H*x(:))./Hs;
xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
% xPhys(voidSet) = 0;
xold1 = x(:);
xold2 = x(:);
%% uniform initialization
k = ones(nely*nelx, 1)*1.0e-2;
%% gradient
g0 = 1.0/nely;
%%
k0 = k;
figure
colormap(parula); imagesc(reshape(k(:)/max(k(:)), nely, nelx)); axis equal; axis tight; axis off; drawnow;
% k(voidSet) = 0;
malpha = 0.1/((nely)^2);
xold1 = [xold1; k(:); g0];
xold2 = [xold2; k(:); g0];
mcon = [];
% change = 1.0;
% erro  = 1.0;
% var_change = [];

while loop < nloop 
    loop = loop+1;
    
    if  mod(loop, 30) == 0 && rou < 100
        rou = rou + 10;
    end
    
    if mod(loop, 20) == 0 && beta < 100
        beta = beta + 2;
    end
    
     %% smoothness
    mstart = [nely+1:nely+1:(nely+1)*(nelx+1)];
    mend = [];
    v = [-1 -1; 1 -1; 1 1; -1 1];
    xxx = ones(nely*nelx, 1);
    rk = k(:);
    [KET]=PreKTQ4_new(v);
    KET1 = speye(size(KET));
    sKT = reshape(KET(:)*rk(:)',16*nelx*nely,1);
    KT = sparse(iKT,jKT,sKT);
    K1 = (KT+KT')/2;
    sKT = reshape(KET1(:)*ones(1, nely*nelx),16*nelx*nely,1);
    K2 = sparse(iKT,jKT,sKT);
    A = K1 + malpha*K2;
    aa = [mstart, mend];
    b = sparse(mstart, ones(1,length(mstart)), 1.0, (nely+1)*(nelx+1), 1);
    B = sparse(aa, aa, ones(1,length(aa)), size(A, 1), size(A, 2));
    M = speye(size(A, 1), size(A, 2));
    M(aa, aa) = 0;
    TA = (B + M*A)';
    T = (B+M*A)\b;
    v = TT'*T;
    tPhys = 1 - v;
    
    %% sensitivity
    II = edofMatFT';
    II = II(:);
    JJ = repmat([1:size(k, 1)]', 1, 4)';
    JJ = JJ(:);
    
    %  %% for density
    % S = ((T(edofMatFT)*KET).*repmat(k(:), 1, 4))';
    % dAdx = M*sparse(II, JJ, S(:));

    %% for conductivity
    S = ((T(edofMatFT)*KET))';
    dAdk =  M*sparse(II, JJ, S(:));
   
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    dv = ones(nely,nelx);
    dx = beta * (1-tanh(beta*(xTilde-eta)).*tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    dv(:) = H*(dv(:).*dx(:)./Hs);
    
    %% objective 
    [c, dc] = Cal_c_ce_whole(nelx, nely, KE, xPhys, Emin, Emax, penal, freedofs, ff);
    f0val = c;
    dc = dc(:);
    dc(:) = H*(dc(:).*dx(:)./Hs);
    df0dx = [dc(:); zeros(nely*nelx+1, 1)]; 

    %% volume constraint
    fval = sum(sum(xPhys)) / (tolvol*volfrac) - 1;
    dfdx = [dv(:)' / (tolvol*volfrac), zeros(1, nely*nelx+1)]; 

    %% layer thickness constraint
    % D = diag(xPhys(:));
    % [G, dGt] = compute_gradient_element(T, F);
    % A = D*LM*G;
    % sk = 1000;
    % B = A.^2;
    % obj = sum(B);
    % fval = [fval; sk*(sum(B)-1.0e-6)];
    % dft = full(2*(D*LM*dGt)'*A);
    % alpha = -TA\dft;
    % dft = sk*alpha'*dAdk;
    % % dfx = alpha'*dAdx;
    % % dfx = H*(dfx(:).*dx(:)./Hs);
    % dfx = sk*full(2*(LM*G).*A);
    % dfx = H*(dfx(:).*dx(:)./Hs);
    % dfdx = [dfdx; dfx(:)', dft(:)'];

    %%
    % D = diag(xPhys(:));
    % [G, dGt] = compute_gradient_element(T, F);
    % sk = 0.02;
    % fval = [fval; sk*sum((D*(G-g0)).^2)];
    % 
    % dft = full(2*(D*dGt)'*(G-g0));
    % alpha = -TA\dft;
    % dft = sk*alpha'*dAdk;
    % 
    % dfx = sk*full(2*(D*(G-g0)).*(G-g0));
    % dfx = H*(dfx(:).*dx(:)./Hs);
    % 
    % dfdx = [dfdx; dfx(:)', dft(:)', -2*sk*sum(D*(D*(G-g0)))];

      
    [G, dGt] = compute_gradient_element(T, F);
    sk = 0.05;
    fval = [fval; sk*sum((G-g0).^2)];

    dft = full(2*(dGt)'*(G-g0));
    alpha = -TA\dft;
    dft = sk*alpha'*dAdk;

    dfdx = [dfdx; zeros(1, nely*nelx), dft(:)', -2*sk*sum(G-g0)];

    %% UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    move = 0.02;
    kmove = 0.01;
    xmin=max(0.0, x(:)-move);
    xmax=min(1, x(:)+move);
    gmove = 0.05;
    xmin = [xmin; max(k-kmove, 0); max(g0-gmove, 0.0)];
    xmax = [xmax; min(k+kmove, 1.0); min(g0+gmove, 1.0)];
    xval = [x(:); k(:); g0];

    %%
    percent = 1 / nStage;
    tP = linspace(0, 1, nStage+1);
    for i = 1 : nStage - 1
        %%
        ti = tP(i+1);
        ft = 1 - (tanh(rou*ti) + tanh(rou*(tPhys - ti)))/(tanh(rou*ti) + tanh(rou*(1-ti)));
        dfdt = -(rou*(tanh(rou*(tPhys - ti)).^2 - 1))/(tanh(rou*(ti - 1)) - tanh(rou*ti));

        xtJoint = xPhys(:).*ft(:);
        fval = [fval; sum(xtJoint(:))/(tolvol*volfrac) - i*percent];

        % for time
        W=xPhys(:).*dfdt(:);
        ss = TT*W;
        alpha = TA\ss;
        dft = alpha'*dAdk/(tolvol*volfrac);

        % for rho
        dfx1 = H*(ft(:).*dx(:)./Hs);
        % dfx2 = alpha'*dAdx;
        % dfx2 = H*(dfx2(:).*dx(:)./Hs);
        dfx = (dfx1)/(tolvol*volfrac);

        %
        dfdx = [dfdx; dfx(:)', dft(:)', 0];

        %
        fval = [fval; -sum(xtJoint(:))/(tolvol*volfrac) + i*percent - 1.0e-2];
        dfdx = [dfdx; -dfx(:)', -dft(:)', 0];
    end

    mcon = [mcon; [f0val; fval]'];
    m=length(fval);
    mdof = 1:m;
    n = size(dfdx, 2);
    
    a0 = 1;
    a = zeros(m,1);
    c_ = ones(m,1)*1000;
    d = zeros(m,1);
    
    [xmma, ymma, zmma, lam, xsi, eta_, mu, zet, s, low, upp] = ...
        mmasub(m, n, loop, xval, xmin, xmax, xold1, xold2,...
        f0val, df0dx, fval(mdof), dfdx(mdof,:),low, upp, a0, a, c_, d);
    
    xold2 = xold1;
    xold1 = xval;

    % changex = max(abs(xmma(1:nely*nelx)-xval(1:nely*nelx)));
    % changek = max(abs(xmma(nely*nelx+1:end)-xval(nely*nelx+1:end)));
    % var_change = [var_change; changex changek];
    x = xmma(1:nely*nelx);
    x(voidSet) = 0;
    xTilde(:) = (H*x(:))./Hs;
    xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    % xPhys(voidSet) = 0; 
    
    %%
    k = xmma(nely*nelx+1 : end-1);
    % k(voidSet) = 0;
    g0 = xmma(end);

    %% PRINT RESULTS
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%4.4f',f0val) ...
        ' Vol.: ' sprintf('%4.3f',sum(sum(xPhys))/(tolvol*volfrac)) ...
        ' Cons.: ' sprintf('%6.3f',fval(2:end))]);

%     if mod(loop, 10) == 0
%         figure(1);
%         colormap(gray); imagesc(-reshape(xPhys(:), nely, nelx), [-1 0]); axis equal; axis tight; axis off; drawnow;
%         hold on
% %         set(gcf,'outerposition',get(0,'screensize'));
% 
% %         filename = sprintf(['density' '\\iter%i.png'], loop);
% %         export_fig(filename);
% 
%         figure(2);
%         colormap(parula); imagesc(reshape(tPhys(:), nely, nelx)); axis equal; axis tight; axis off; drawnow;
% %         set(gcf,'outerposition',get(0,'screensize'));
% 
% %         filename = sprintf(['time' '\\iter%i.png'], loop);
% %         export_fig(filename);
%     end
end
