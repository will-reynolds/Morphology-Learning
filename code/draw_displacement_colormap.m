% function draw_displacement(xPhys, 
close all
load U1.mat
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
figure(1);
% axis([-0.5 nelx+1 0 nely+5.5]); % for planar
% % axis([-0 nelx+1 0 nely+1]); % for optimized
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

dis1 = V1 - V;
dis1 = sqrt(sum(dis1.^2, 2));
% dis1 = dis1 / max(dis1);

%
mMove = reshape(UU2, 2, [])';
V2 = V + mMove;
dis2 = V2 - V;
dis2 = sqrt(sum(dis2.^2, 2));
% dis2 = dis2 / max(dis2);

dis0 = unique([dis1; dis2]);
dis = dis0 / max(dis0);

color = jet(100);
mycolor = zeros(length(dis1), 3);

%
for i = 1 : length(dis1)
    [xx1, yy] = find(dis0 == dis1(i));
    
    a = floor(dis(xx1)*100);
    if a == 0
        a = 1;
    end
    mycolor(i, :) = color(a, :);
end

% patch(mX(M(:,xx)),mY(M(:,xx)),'r', 'FaceAlpha', 0.3);
% figure(1)
patch('Faces', M(:, xx)', 'Vertices', V1, 'FaceVertexCData', mycolor, 'Facecolor', 'interp');
axis equal;
% axis([-0.5 nelx+1 0 nely+5.5]); % for planar
axis([-0 nelx+1 0 nely+1]); % for optimizedaxis off
cameratoolbar;

%
mycolor = zeros(length(dis2), 3);
figure(2)
hold on
patch(X(M(:,xx)),Y(M(:,xx)),[0.8, 0.8, 0.8],'FaceAlpha', 0.3);
[xx1, yy1] = find(Mask(:) == 1);

for i = 1 : length(dis2)
    [xx1, yy] = find(dis0 == dis2(i));
    
    a = floor(dis(xx1)*100);
    if a == 0
        a = 1;
    end
    mycolor(i, :) = color(a, :);
end

% patch(mX(M(:,xx)),mY(M(:,xx)),'r', 'FaceAlpha', 0.3);

patch('Faces', M(:, xx)', 'Vertices', V2, 'FaceVertexCData', mycolor, 'Facecolor', 'interp');
axis equal;
axis([-0.5 nelx+1 0 nely+5.5]); % for planar
% axis([-0 nelx+1 0 nely+1]); % for optimizedaxis off
cameratoolbar;

% mMove(1, :)
% mMove(nelx*(nely+1)+1, :)
% val = norm(mMove(1, :) - mMove(nelx*(nely+1)+1, :))^2

% UU'*Q*UU
% % 
% figure
% ss = 1;
% x=ss:nloop;
% plot(x, mcon(ss:end,1), 'r-','linewidth',4)
% set(gca, 'fontsize', 30)

% end