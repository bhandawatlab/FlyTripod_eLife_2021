function Rn=getLegLengthEditor(data)

%data = data.frames;


FigH = figure('position',[100 100 size(data,2)+100 size(data,1)+600]);
% axes('XLim', [0 4*pi], 'units','pixels', ...
%      'position',[100 50 200 200], 'NextPlot', 'add');
% x     = linspace(0, 4*pi, 400);
% y     = sin(x);
dcm_obj = datacursormode(FigH);
set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on')

first = 1;
imgH = imshow(data(:,:,:,first));
TextH = uicontrol('style','text',...
    'position',[100 40 40 15]);
TextH.String='1';
SliderH = uicontrol('style','slider','position',[1 10 1200 20],...
    'min', 1, 'max', size(data,4),'Value',1);
addlistener(SliderH, 'Value', 'PostSet', @callbackfn);
movegui(FigH, 'center')


Rnat=zeros(5,1);
for j=1:5
i=1;
joints=zeros(8,2);
while true
    
    while true
        w = waitforbuttonpress;
        if w == 1
            num = get(gcf, 'CurrentCharacter');
            num = str2double(num)
            if ismember(num,[0 1])
                break
            end
        end
    end
    
    if num == 1
        disp(['registered: ' num2str(i)])
        c_info = getCursorInfo(dcm_obj);
        position = c_info.Position;
    else
        disp('exiting prematurely')
        break
    end
    
    joints(i,:)=position;
    
    i=i+1;
    
    if i == 9
        i=1;
        break
    end
end


%joins row 1~4 bottom view. row 5~8 top view
%From tip to body
%
joints3D=zeros(4,3);

joints3D(1,1)=joints(1,1);
joints3D(2,1)=joints(2,1);
joints3D(3,1)=joints(3,1);
joints3D(4,1)=mean(joints([4 8],1));

joints3D(1,2)=joints(1,2);
joints3D(2,2)=joints(2,2);
joints3D(3,2)=joints(3,2);
joints3D(4,2)=joints(4,2);

joints3D(1,3)=joints(5,2);
joints3D(2,3)=joints(6,2);
joints3D(3,3)=joints(7,2);
joints3D(4,3)=joints(8,2);
joints3D(:,3)=joints3D(1,3)-joints3D(:,3);

% figure
% plot3(joints3D(:,1),joints3D(:,2),joints3D(:,3))
% xlabel('x')
% ylabel('y')
% zlabel('z')

segmentLengths = sqrt(sum(diff(joints3D,1,1).^2))*24/1984;
Rnat(j) = sum(segmentLengths);
disp('Move to next frame')
disp(['currently ' num2str(j) '/5'])
end
Rn = zeros(1,2);
Rn(1)=mean(Rnat);
Rn(2)=std(Rnat);

    function callbackfn(source, eventdata)
        num          = round(get(eventdata.AffectedObject, 'Value'));
        imgH.CData = data(:,:,:,first*num);
        
        TextH.String = num2str(num);
    end

end