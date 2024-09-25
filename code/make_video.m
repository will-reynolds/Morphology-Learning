[file, pathname] = uigetfile('*.png', 'multiselect', 'on');
filename = './movie/curved_2';
% 
aviobj = VideoWriter([filename '.mp4'], 'MPEG-4');
aviobj.FrameRate = 1;
aviobj.Quality = 100;
open(aviobj)
t = 0;
for i = 1 : length(file)
    
%         frame = imread([pathname, filename{i}]);
    im=imread([pathname, file{i}]);
%     if i == 1
%         t = size(im);
%     end
%     
%     if i >= 106
%         a = zeros(t);
%         b = size(im);
%         a(1:b(1), 1:b(2), 1:b(3)) = im;
%         im = a/255;
%     end
%     [I,map]=rgb2ind(im,65535);
%     I= uint8(I);
%     if i==1
%         imwrite(I,map,[filename '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);%FIRST
%     else
%         imwrite(I,map,[filename '.gif'],'gif','WriteMode','append','DelayTime',0.1);
%     end
    writeVideo(aviobj, im);
end

close(aviobj);