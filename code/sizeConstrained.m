function [con, dft, dfx] = sizeConstrained(xPhys, tPhys, nelx, nely, nStage, rou, p)
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
V = [X(:), Y(:)];
F = M';
VE = (V(F(:, 1), :) + V(F(:, 2), :) + V(F(:, 3), :) + V(F(:, 4), :))/4;

%%
tP = linspace(0, 1, nStage+1);
dfx = [];
dft = [];
con = [];
for i = 1 : nStage
    %%
    ti = tP(i+1);
    tj = tP(i);
    ft = 1 - (tanh(rou*ti) + tanh(rou*(tPhys - ti)))/(tanh(rou*ti) + tanh(rou*(1-ti)));
    ft_prev = 1 - (tanh(rou*tj) + tanh(rou*(tPhys - tj)))/(tanh(rou*tj) + tanh(rou*(1-tj)));

    dfdt = -(rou*(tanh(rou*(tPhys - ti)).^2 - 1))/(tanh(rou*(ti - 1)) - tanh(rou*ti));
    dfdt_prev = -(rou*(tanh(rou*(tPhys - tj)).^2 - 1))/(tanh(rou*(tj - 1)) - tanh(rou*tj));

    if i == 1
        mask = repmat(xPhys(:) .* ft(:), 1, 2) .* VE;
        w = sum(xPhys(:) .* ft(:));
        cen = sum(mask, 1) / w;
        var = (VE - repmat(cen, size(VE, 1), 1)) ;
        dis = sum(var.*var, 2).*xPhys(:) .* ft(:);
        con = [con; sum(dis.^p)^(1/p) - nely/2];

        %
        a1 = sum(dis.^p)^(1/p-1);
        a2 = sum(var.*var, 2).*xPhys(:).*dfdt(:).*dis.^(p-1);
        vec = sum(repmat(xPhys(:).*ft(:).*dis.^(p-1), 1, 2).*var, 1);
        a3 = sum((w*VE-repmat(sum(mask, 1), size(VE, 1), 1)).*repmat(vec, size(VE, 1), 1), 2);
        a4 = -2/w^2*xPhys(:).*dfdt(:).*a3;
        dft = [dft; [a1.*(a2 + a4)]'];
        %
        a2 = sum(var.*var, 2).*ft(:).*dis.^(p-1);
        a4 = -2/w^2*ft(:).*a3;
        dfx = [dfx; [a1.*(a2 + a4)]'];
    else 
        mask = repmat(xPhys(:) .* (ft(:)-ft_prev(:)), 1, 2) .* VE;
        cen = sum(mask, 1) / sum(xPhys(:) .* (ft(:)-ft_prev(:)));
        var = (VE - repmat(cen, size(VE, 1), 1)) ;
        dis = sum(var.*var, 2).*xPhys(:) .* (ft(:)-ft_prev(:));
        con = [con; sum(dis.^p)^(1/p) - nely/2];
        %
        a1 = sum(dis.^p)^(1/p-1);
        a2 = sum(var.*var, 2).*xPhys(:).*(dfdt(:)-dfdt_prev(:)).*dis.^(p-1);
        vec = sum(repmat(xPhys(:).*(ft(:)-ft_prev(:)).*dis.^(p-1), 1, 2).*var, 1);
        a3 = sum((w*VE-repmat(sum(mask, 1), size(VE, 1), 1)).*repmat(vec, size(VE, 1), 1), 2);
        a4 = -2/w^2*xPhys(:).*(dfdt(:)-dfdt_prev(:)).*a3;
        dft = [dft; [a1.*(a2 + a4)]'];
        %
        a2 = sum(var.*var, 2).*(ft(:)-ft_prev(:)).*dis.^(p-1);
        a4 = -2/w^2*ft(:).*a3;
        dfx = [dfx; [a1.*(a2 + a4)]'];
    end

end

end