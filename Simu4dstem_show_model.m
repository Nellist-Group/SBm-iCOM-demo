function  Simu4dstem_show_model( CrysData, Nf )
%SHOWMODEL Show Model, less than 7 kinds of elements in the model
%   CrysData = �������룬ԭ������xyz��꣬��λA
%   Nf = ���Ϊfigure(Nf)

%ShowModelVersion = 170727

CrysData = sortrows(CrysData);      %��CrysData����ԭ����������
Elements = unique(CrysData(:,1))';  %ȷ��ģ�����м���ԭ��
Elements = [Elements; zeros(size(Elements))];   
%Elements��һ��Ϊģ�����е�ԭ������ڶ���Ϊ���ԭ����ģ���е�����

ColorOrder = 'rgbcmwk';

if size( Elements,2 ) > size(ColorOrder,2)
    error('Function ShowModel error: need more ColorOrder');
end

figure(Nf);
clear gcf;
hold on;
axis equal;
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
view(0,90);

count = 1;
for ii = 1 : size(Elements,2)       
    %ѭ������ÿһ��ԭ����ģ���е����������Ҳ��ģ�͵�����ԭ��
    Elements(2,ii) = sum(CrysData(:,1) == Elements(1,ii));
    ElmPlot = CrysData( count: count + Elements(2,ii) - 1,:);
    count = count + Elements(2,ii);
    plot3(ElmPlot(:,2),ElmPlot(:,3),ElmPlot(:,4),[ColorOrder(ii),'.']);
end

% load('ElementData.mat');
% figure(Nf + 1);
% %axes(haxes);    %ѡ����ͼaxes
% hold on;
% axis equal;
% axis off;
% light('Position',[1 1 2]);  % ����
% light('Position',[-1,-1,-2]);
% [x,y,z] = sphere;   % ��ɱ�׼������
% for ii = 1:size(CrysData,1)    %��ÿһ�У���ÿһ��ԭ�ӻ�ͼ
%     r = ElementData(CrysData(ii,1),3);
%     color = ElementData(CrysData(ii,1),4:6);
%     surface(CrysData(ii,2)+r*x,CrysData(ii,3)+r*y,...
%         CrysData(ii,4)+r*z,...
%         'FaceColor',color,'EdgeColor','none','FaceLighting','gouraud');
% end

end

