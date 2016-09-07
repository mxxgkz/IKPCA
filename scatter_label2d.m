function scatter_label2d(Z,title_text,dd,tag,psize,color)
% Plot scatter plots and label each data (coded with color which varied in 2-dimension)
% Input:
% Z: matrix of input data points with rows as points (points are
% 2-dimensional
% title_text: The title text used in plot
% dd: the factor of shifted translation for labels
% tag: A tag to decide whether indices of data points should be marked
% besides those data points (0: no indices; 1: indices)
% color: color map for each data points (rows of Z)

N = size(Z,1); %Number of data points

figure();
max_range=max(range(Z,1));
axis equal;
scatter(Z(:,1),Z(:,2),psize,color,'filled'); %Scatter plot
axis equal;
range_Z = range(Z,1); %Range of each coordinates
a = (1:N)'; b = num2str(a); c = cellstr(b); %Labels
if tag == 1
    text(Z(:,1)+dd*range_Z(1),Z(:,2)+dd*range_Z(2),c); %Label those data points    
end
title(title_text);
end