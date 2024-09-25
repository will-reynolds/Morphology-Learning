[file, pathname] = uigetfile('*.png', 'multiselect', 'on');
filename = './movie/robot';

for i = 1 : length(file)
    
    im=imread([pathname, file{i}]);
%     if size(im ,3) == 1 
%         im = repmat(im ,1,1,3);
%     end
    [I,map]=rgb2ind(im,256);
    I= uint8(I);
    if i==1
        imwrite(I,map,[filename '.gif'],'gif', 'Loopcount',inf,'DelayTime', 1);%FIRST
    else
        imwrite(I,map,[filename '.gif'],'gif','WriteMode','append','DelayTime', 1);
    end

end
