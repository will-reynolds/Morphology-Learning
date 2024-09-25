% function draw_displacement(xPhys, 
close all
% U=mU{i};
% Mask=mTiming{i};
Mask = xPhys;
%% element and point relationship
M = zeros(4, nely*nelx);
for x = 1 : nelx
    for y = 1 : nely
        i = nely*(x - 1) + y;
        M(:, i) = [(nely+1)*x + y; (nely+1)*(x-1) + y; (nely+1)*(x-1)+y + 1; (nely+1)*x + y + 1];
    end
end

%% draw
X = repmat(0.5:nelx+.5, nely+1, 1);
X = flipud(X);
Y = repmat([0.5:nely+.5]', 1, nelx+1);
Y = flipud(Y);
V = [X(:), Y(:)];
F = M';

F_c = (V(F(:, 1), :) + V(F(:, 2), :) + V(F(:, 3), :) + V(F(:, 4), :)) / 4;

[xx, yy] = find(xPhys(:) > 0.1);
% xx = 1 : nely*nelx;
%% draw
figure;
axis([0 nelx+1 0 nely + 1]);
hold on
patch(X(M(:,xx)),Y(M(:,xx)),[0.8, 0.8, 0.8],'FaceAlpha', 0.3);
[xx1, yy1] = find(Mask(:) == 1);
hold on
% patch(X(M(:,xx1)),Y(M(:,xx1)),[0., 0., 0.],'FaceAlpha', 0.6);
%     patch('Faces', F(1:a, :), 'Vertices', V, 'FaceCoaxlor', [0.8, 0.8, 0.8]);
hold on
%% deformed
mMove = reshape(UU, 2, [])';
V1 = V + mMove;
mX = V1(:, 1);
mX = reshape(mX, nely+1, nelx+1);
mY = V1(:, 2);
mY = reshape(mY, nely+1, nelx+1);

patch(mX(M(:,xx)),mY(M(:,xx)),'r', 'FaceAlpha', 0.3);

axis equal;
axis([-20 nelx+2 -20 nely + 8]);
axis off
cameratoolbar;
% mMove(1, :)
% mMove(nelx*(nely+1)+1, :)
% val = norm(mMove(1, :) - mMove(nelx*(nely+1)+1, :))^2
hold on 
% draw_points(V1(nn, :), 30, [0 1 1]);
hold on 
% draw_points(V1([1 ], :), 120, [0 1 1]);
UU'*Q*UU/length(nn1)
% % 
set(gcf, 'color', 'w');
figure
ss = 1;
x=ss:nloop;
plot(x, mcon(ss:end,1), 'r-','linewidth',4)
set(gca, 'fontsize', 30, 'linewidth', 3)
%%

% end