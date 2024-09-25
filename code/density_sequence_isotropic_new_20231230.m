clear
close all
tic
% A = imread('80.png');
% A = rgb2gray(A)/255;
% A = A > 0.1;
% xPhys= 1-A;
% [xx, yy] = find(xPhys(:) == 0);
% voidSet = xx;
nelx = 210;
nely = 140;
xPhys = ones(nely, nelx);
voidSet = [];
nStage = 20;
nloop = 500;
volfrac = 0.5;
% [nely, nelx] = size(xPhys);
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
Y = flipud(Y);
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
rmin = 5;     % density filter radiu

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
% H(:, [fixedSet]) = 0;
Hs = sum(H,2);

% %% anisotropic
% D = -[0.5 0 0.5; 0 0.5 0.5; -0.5 0 0.5; 0 0.5 -0.5; -0.5 0 -0.5; 0 -0.5 -0.5; 0.5 0 -0.5; 0 -0.5 0.5]';
% S = Emax*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]/(1-nu^2);
% I = zeros(4*nelx*nely, 1); 
% J = zeros(4*nelx*nely, 1); 
% Sx = zeros(4*nelx*nely, 1); 
% Sy = zeros(4*nelx*nely, 1);
% inelastic_str = [-0.1, 0, 0]';
% strain0 = D'*S*inelastic_str;
% % f_val = inelastic_str*Emax /(2*(1 - nu));
% k = 0;
% for x = 1 : nely
%     for y = 1 : nelx
%         I(k+1:k+4) = [(y-1)*(nely+1)+x+1; y*(nely+1)+x+1; y*(nely+1)+x; (y-1)*(nely+1) + x];
%         J(k+1:k+4) = [(y-1)*nely + x; (y-1)*nely + x; (y-1)*nely + x; (y-1)*nely + x];
%         Sx(k+1:k+4) = strain0(1:2:end);
%         Sy(k+1:k+4) = strain0(2:2:end);
%         k = k + 4;
%     end
% end 
% % Cx = f_val*sparse(I, J, Sx);
% % Cy = f_val*sparse(I, J, Sy);
% Cx = sparse(I, J, Sx);
% Cy = sparse(I, J, Sy);
% CC = sparse(2*size(Cx, 1), size(Cx, 2));
% CC(1 : 2 : end, :) = Cx;
% CC(2 : 2 : end, :) = Cy;
% CC = sparse(CC);

%% isotropy
%%
I = zeros(4*nelx*nely, 1); 
J = zeros(4*nelx*nely, 1); 
Sx = zeros(4*nelx*nely, 1); 
Sy = zeros(4*nelx*nely, 1);
inelastic_str = -0.1*60/nelx;
f_val = inelastic_str*Emax /(2*(1 - nu));
k = 0;
for x = 1 : nely
    for y = 1 : nelx
        I(k+1:k+4) = [(y-1)*(nely+1) + x; (y-1)*(nely+1)+x+1; y*(nely+1)+x; y*(nely+1)+x+1];
        J(k+1:k+4) = [(y-1)*nely + x; (y-1)*nely + x; (y-1)*nely + x; (y-1)*nely + x];
        Sx(k+1:k+4) = [-1; -1; 1; 1];
        Sy(k+1:k+4) = [1; -1; 1; -1];
        k = k + 4;
    end
end
Cx = f_val*sparse(I, J, Sx);
Cy = f_val*sparse(I, J, Sy);
CC = sparse(2*size(Cx, 1), size(Cx, 2));
CC(1 : 2 : end, :) = Cx;
CC(2 : 2 : end, :) = Cy;
CC = sparse(CC);

%%
nn = [1, (nely+1)*nelx + 1]; 
% nn = 1 : nely+1 : (nely+1)*nelx+1;
Q = sparse(2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1));
for i = 1 : length(nn)
    Q1 = sparse(2*(nely+1)*(nelx+1),2*(nely+1)*(nelx+1));
    Q1(1, nn*2) = -1/length(nn);
    Q1(1, 2*nn(i)) = Q1(1, 2*nn(i)) + 1;
    Q1 = Q1'*Q1;
    Q = Q + Q1;
end
hold on
draw_points(V(nn, :), 20, [1 0 0]);
cameratoolbar

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
%
fixdeNode_thermal = nely+1:nely+1:(nely+1)*(nelx+1);
[xx, yy] = find(xPhys == 1);
aa = [(yy-1)*(nely+1)+xx, (yy-1)*(nely+1)+xx+1, yy*(nely+1)+xx, yy*(nely+1)+xx+1];
aa = unique(aa(:));
bb = 1 : (nely+1)*(nelx+1);
aa = setdiff(bb, aa);
fixedNode = [fixdeNode_thermal, aa];
fixeddofs = unique([2*fixdeNode_thermal, 2*fixdeNode_thermal-1]);
alldofs = 1 : 2*(nely+1)*(nelx+1);
freedofs_thermal = setdiff(alldofs, fixeddofs);

%%
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

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
rou = 10;
loop = 0;
beta = 1;       % beta continuation
eta = 0.5;      % projection threshold, fixed at 0.5

x = repmat(volfrac,nely,nelx);
xTilde = x;
xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
xPhys(voidSet) = 0;
xold1 = x(:);
xold2 = x(:);
% k = ones(nely*nelx, 1)*1.0e-3;
% k = rand(nely*nelx, 1)*1.0e-3;
k = ones(nely*nelx, 1)+10^3;
k0 = k;
figure
colormap(parula); imagesc(reshape(k(:)/max(k(:)), nely, nelx)); axis equal; axis tight; axis off; drawnow;
k(voidSet) = 0;
% malpha = 0.1/((2*nelx)^2);
xold1 = [xold1; k(:)];
xold2 = [xold2; k(:)];
mcon = [];

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
    v = [-1 -1; 1 -1; 1 1; -1 1];
    rk = xPhys(:).*k(:);
    [iKT,jKT,KET,TT,edofMatFT]=PreKTQ4(nelx, nely, v);
    TT = sparse(TT);
    sKT = reshape(KET(:)*rk(:)',16*nelx*nely,1);
    KT = sparse(iKT,jKT,sKT);
    KT = (KT+KT')/2;
    xx = 1:size(KT, 1);
    xx = setdiff(xx, mstart);
    yy = xx;
    ss = ones(1,length(xx));
    W = sparse(xx, yy, ss, (nely+1)*(nelx+1), (nely+1)*(nelx+1));
    A = speye(size(KT)) + W*KT;
    TA = A';
    b = zeros((nely+1)*(nelx+1), 1);
    b(mstart) = 1;
    T = A\b;
    v = TT'*T;
    tPhys = 1 - v;
%     tPhys(:) = (H*tPhys(:))./Hs;

%     if loop == 400
%         colormap(summer); imagesc(reshape(1-tPhys(:), nely, nelx)); axis equal; axis tight; axis off; drawnow;
%         figure
%         colormap(parula); imagesc(reshape(tPhys(:), nely, nelx)); axis equal; axis tight; axis off; drawnow;
%     end
    
    %% sensitivity
    II = edofMatFT';
    II = II(:);
    JJ = repmat([1:size(k, 1)]', 1, 4)';
    JJ = JJ(:);
    
     %% for density
    S = ((T(edofMatFT)*KET).*repmat(k(:), 1, 4))';
    dAdx = W*sparse(II, JJ, S(:));

    %% for conductivity
    S = ((T(edofMatFT)*KET).*repmat(xPhys(:), 1, 4))';
    dAdk =  W*sparse(II, JJ, S(:));
   
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    dv = ones(nely,nelx);
    dx = beta * (1-tanh(beta*(xTilde-eta)).*tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    dv(:) = H*(dv(:).*dx(:)./Hs);
%     dv(:) = H*(dv(:)./Hs);
    % dv([fixedSet]) = [];
    
    %% objective 
    [c, dc] = Cal_c_ce_whole(nelx, nely, KE, xPhys, Emin, Emax, penal, freedofs, ff);
    f0val = c;
    dc = dc(:);
    dc(:) = H*(dc(:).*dx(:)./Hs);
%     dc(:) = H*(dc(:)./Hs);
    % dc([fixedSet]) = [];
    df0dx = [dc(:); zeros(nely*nelx, 1)]; 
        
    %% volume constraint
    fval = sum(sum(xPhys))/(tolvol*volfrac) - 1;
    dfdx = [dv(:)'/(tolvol*volfrac), zeros(1, nely*nelx)]; 

%      %% distortion constraint
%     tP = linspace(0, 1, nStage + 1);
%     UU = zeros(2*(nely+1)*(nelx+1), 1);
%     dKdtU = cell(nStage, 1);
%     dKdxU = cell(nStage, 1);
%     LL = cell(nStage, 1);
% 
%     for i = 1 : nStage
%         ti = tP(i+1);
%         tj = tP(i);
%         gt = 1-(tanh(rou*tj) + tanh(rou*(tPhys(:) - tj)))/(tanh(rou*tj) + tanh(rou*(1-tj)));
%         dgdt = -(rou*(tanh(rou*(tPhys(:) - tj)).^2 - 1))/(tanh(rou*(tj - 1)) - tanh(rou*tj));
% 
%         ft = 1 - (tanh(rou*ti) + tanh(rou*(tPhys(:) - ti)))/(tanh(rou*ti) + tanh(rou*(1-ti)));
%         dfdt = -(rou*(tanh(rou*(tPhys(:) - ti)).^2 - 1))/(tanh(rou*(ti - 1)) - tanh(rou*ti));
% 
%         if i == nStage
%             xtJoint = xPhys(:);
%         else
%             xtJoint = xPhys(:).*ft(:);
%         end
% 
%         %% for force
%         % xtfJoint = xPhys.*ft.*(1-gt);
%         if i == 1
%             xtfJoint = xPhys(:).*ft(:);
%         elseif i == nStage
%             xtfJoint = xPhys(:).*(1-gt(:));
%         else
%             xtfJoint = xPhys(:).*(ft(:)-gt(:));
%         end
% 
%         h = 3;
%         xt = xtfJoint(:).^h;
%         FF = CC*xt(:);
%         FF = sparse(FF);
% 
%         %
%         sK = reshape(KE(:)*(E0+xtJoint(:)'.^penal*(Emax-Emin)),64*nely*nelx,1);
%         K = sparse(iK,jK,sK);
%         K = (K+K')/2;
%         LL{i} = K;
% 
%         U = zeros(2*(nely+1)*(nelx+1), 1);
%         U(freedofs_thermal) = K(freedofs_thermal, freedofs_thermal)\FF(freedofs_thermal);
% 
%         UU = UU + U;
%         a = KE*U(edofMat');
% 
%         clear K U
% 
%         if i == 1
%             dKdtU{i, 1} = edofMat;
%             dKdtU{i, 2} = repmat([1:size(a, 2)]', 1, size(edofMat, 2));
%             dKdxU{i, 1} = edofMat;
%             dKdxU{i, 2} = repmat([1:size(a, 2)]', 1, size(edofMat, 2));
%         else
%             dKdtU{i, 1} = [];
%             dKdtU{i, 2} = [];
%             dKdxU{i, 1} = [];
%             dKdxU{i, 2} = [];
%         end
% 
%         if i == nStage
%             dKdtU{i, 3} = a'*0;
%             b = penal*xtJoint(:).^(penal-1)*(Emax-Emin);
%             dKdxU{i, 3} =  (a*sparse(1:length(b), 1:length(b), b))';
%         else
%             b = penal*xPhys(:).*dfdt(:).*xtJoint(:).^(penal-1)*(Emax-Emin);
%             dKdtU{i, 3} =  (a*sparse(1:length(b), 1:length(b), b))';
% 
%             b = penal*ft(:).*xtJoint(:).^(penal-1)*(Emax-Emin);
%             dKdxU{i, 3} =  (a*sparse(1:length(b), 1:length(b), b))';
% 
%         end
% 
%         clear b a
%     end
% 
%     epsilon = 1.0e-3;
%     fval = [fval; UU'*Q*UU/length(nn)- epsilon];
% 
%     %% sensitivity
%     dct = zeros(1, nely*nelx);
%     dcx = zeros(1, nely*nelx);
%     for i = 1 : nStage
%         ti = tP(i+1);
%         tj = tP(i);
%         gt = 1-(tanh(rou*tj) + tanh(rou*(tPhys - tj)))/(tanh(rou*tj) + tanh(rou*(1-tj)));
%         dgdt = -(rou*(tanh(rou*(tPhys - tj)).^2 - 1))/(tanh(rou*(tj - 1)) - tanh(rou*tj));
% 
%         ft = 1 - (tanh(rou*ti) + tanh(rou*(tPhys - ti)))/(tanh(rou*ti) + tanh(rou*(1-ti)));
%         dfdt = -(rou*(tanh(rou*(tPhys - ti)).^2 - 1))/(tanh(rou*(ti - 1)) - tanh(rou*ti));
% 
%         if i == nStage
%             xtJoint = xPhys(:);
%         else
%             xtJoint = xPhys(:).*ft(:);
%         end
% 
%         %% for force
%         if i == 1
%             xtfJoint = xPhys(:).*ft(:);
%         elseif i == nStage
%             xtfJoint = xPhys(:).*(1-gt(:));
%         else
%             xtfJoint = xPhys(:).*(ft(:)-gt(:));
%         end
% 
%         L = LL{i};
%         FF = -2*Q*UU/length(nn);
% %         FF = sparse(FF);
% 
%         U = zeros(2*(nely+1)*(nelx+1), 1);
%         U(freedofs_thermal) = L(freedofs_thermal, freedofs_thermal)\FF(freedofs_thermal);
% %         U(freedofs) = L'\ (L\FF(freedofs));
%         II = dKdtU{1, 1}; II = II(:);
%         JJ = dKdtU{1, 2}; JJ = JJ(:);
%         SS = dKdtU{i, 3}; SS = SS(:);
%         Kt = sparse(II, JJ, SS);
%         II = dKdxU{1, 1}; II = II(:);
%         JJ = dKdxU{1, 2}; JJ = JJ(:);
%         SS = dKdxU{i, 3}; SS = SS(:);
%         Kx = sparse(II, JJ, SS);
% 
%         %
%         if i == 1
%             dFt = CC(freedofs_thermal, :) * sparse(1:length(xtfJoint(:)), 1:length(xtfJoint(:)), (h*xtfJoint(:).^(h-1).*xPhys(:).*(dfdt(:)))');
%             dFx = CC(freedofs_thermal, :) * sparse(1:length(xtfJoint(:)), 1:length(xtfJoint(:)), (h*xtfJoint(:).^(h-1).*ft(:))');
%         elseif i == nStage
%             dFt = CC(freedofs_thermal, :) * sparse(1:length(xtfJoint(:)), 1:length(xtfJoint(:)), (h*xtfJoint(:).^(h-1).*xPhys(:).*(-dgdt(:)))');
%             dFx = CC(freedofs_thermal, :) * sparse(1:length(xtfJoint(:)), 1:length(xtfJoint(:)), (h*xtfJoint(:).^(h-1).*(1-gt(:)))');
%         else
%             dFt = CC(freedofs_thermal, :) * sparse(1:length(xtfJoint(:)), 1:length(xtfJoint(:)), (h*xtfJoint(:).^(h-1).*xPhys(:).*(dfdt(:)-dgdt(:)))');
%             dFx = CC(freedofs_thermal, :) * sparse(1:length(xtfJoint(:)), 1:length(xtfJoint(:)), (h*xtfJoint(:).^(h-1).*(ft(:)-gt(:)))');
%         end
% 
%         dct = dct + U(freedofs_thermal)'*dFt - U(freedofs_thermal)'*Kt(freedofs_thermal, :);
%         dcx = dcx + U(freedofs_thermal)'*dFx - U(freedofs_thermal)'*Kx(freedofs_thermal, :);
%     end
%     clear dKdtU Kt Kx dFt dFx II JJ SS
% 
%     %%
%      ss = TT*dct';
%     alpha = -TA\ss;
%     dck = alpha'*dAdk;
% 
%     dcx1 = H*(dcx(:).*dx(:)./Hs);
%     dcx2 = alpha'*dAdx;
%     dcx2 = H*(dcx2(:).*dx(:)./Hs);
%     dcx = (-dcx1+dcx2);
%     % dcx([fixedSet]) = [];
% 
%     dfdx = [dfdx; dcx(:)', dck(:)'];
    
    %% UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    move = 0.02;  
    xmin=max(0.0, x(:)-move);
    xmax=min(1, x(:)+move);
    kmove = (10^7-10^2)/50;
    xmin = [xmin; max(k-kmove, 10^2)];
    xmax = [xmax; min(k+kmove, 10^7)];
    xval = [x(:); k(:)];

    %%
    percent = 1 / nStage;
    tP = linspace(0, 1, nStage+1);
    for i = 1 : nStage
        %%
        ti = tP(i+1);
        ft = 1 - (tanh(rou*ti) + tanh(rou*(tPhys - ti)))/(tanh(rou*ti) + tanh(rou*(1-ti)));
        dfdt = -(rou*(tanh(rou*(tPhys - ti)).^2 - 1))/(tanh(rou*(ti - 1)) - tanh(rou*ti));
        
        xtJoint = xPhys(:).*ft(:);
        fval = [fval; sum(xtJoint(:))/(tolvol*volfrac) - i*percent];
   
        % for time
        W = xPhys(:).*dfdt(:);
%         W = H*(W(:)./Hs);
        ss = TT*W;
        alpha = TA\ss;
        dft = alpha'*dAdk/(tolvol*volfrac);
        
        % for rho
        dfx1 = H*(ft(:).*dx(:)./Hs);
        dfx2 = alpha'*dAdx;
        dfx2 = H*(dfx2(:).*dx(:)./Hs);
        dfx = (dfx1+dfx2)/(tolvol*volfrac);
        
        %
        % dfx([fixedSet]) = [];
        dfdx = [dfdx; dfx(:)', dft(:)'];
        
        %
        fval = [fval; -sum(xtJoint(:))/(tolvol*volfrac) + i*percent - 5.0e-2];
        dfdx = [dfdx; -dfx(:)', -dft(:)'];
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
    
    x = xmma(1:nely*nelx);
    x(voidSet) = 0;
    xTilde(:) = (H*x(:))./Hs;
%     xTilde(:) = (H*xTilde(:))./Hs;
    xPhys = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    xPhys(voidSet) = 0; 
    
    %%
    k = xmma(nely*nelx+1 : end);
    k(voidSet) = 0;
    
    %% PRINT RESULTS
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',f0val) ...
        ' Vol.: ' sprintf('%6.3f',sum(sum(xPhys))/(tolvol*volfrac))  ...
        ' Cons.: ' sprintf('%6.3f',fval(2:end))]);
    
    if mod(loop, 10) == 0
        figure(1);
        colormap(gray); imagesc(-reshape(xPhys(:), nely, nelx), [-1 0]); axis equal; axis tight; axis off; drawnow;
        hold on
%         set(gcf,'outerposition',get(0,'screensize'));
        
%         filename = sprintf(['density' '\\iter%i.png'], loop);
%         export_fig(filename);

        figure(2);
        colormap(parula); imagesc(reshape(tPhys(:), nely, nelx)); axis equal; axis tight; axis off; drawnow;
%         set(gcf,'outerposition',get(0,'screensize'));
        
%         filename = sprintf(['time' '\\iter%i.png'], loop);
%         export_fig(filename);
    end
%     
    %     U = UU;
end
toc