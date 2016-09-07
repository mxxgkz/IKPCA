function scatter_label2d(Z,title_text,dd)
% Plot scatter plots and label each data
% Input:
% Z: matrix of input data points with rows as points (points are
% 2-dimensional
% title_text: The title text used in plot
% dd: the factor of shifted translation for labels

N = size(Z,1); %Number of data points
range_factor = 1.2; %Factor multiplied to range of axis of plot
% In case the Z has complex components
Z_real = real(Z);
figure();
hold on;
max_range=range_factor*max(range(Z_real,1)); %max of Range of each coordinates
axis equal;
scatter(Z_real(:,1),Z_real(:,2)); %Scatter plot
axis([mean(Z_real(:,1))-max_range/2,mean(Z_real(:,1))+max_range/2,mean(Z_real(:,2))-max_range/2,mean(Z_real(:,2))+max_range/2]);
a = (1:N)'; b = num2str(a); c = cellstr(b); %Labels
text(Z(:,1)+dd*max_range,Z(:,2)+dd*max_range,c); %Label those data points
title(title_text);
end