function scatter_3d(Z,title_text,color_3D,psize)
% Plot 3-D scatter plot (with equal axis unit and length in the plot)
% Z: input data with columns as different coordinates
% pc1,pc2,pc3: the index of the first three principle components 
% title_text: the title of scatter plot

N = size(Z,1); %Number of data points
range_factor = 1.2; %Factor multiplied to range of axis of plot
% In case the Z has complex components
Z_real = real(Z);
%figure();
hold on;
grid on;
max_range=range_factor*max(range(Z_real,1)); %max of Range of each coordinates
axis equal;
scatter3(Z_real(:,1),Z_real(:,2),Z_real(:,3),psize,color_3D); %Scatter plot
axis([mean(Z_real(:,1))-max_range/2,mean(Z_real(:,1))+max_range/2,...
    mean(Z_real(:,2))-max_range/2,mean(Z_real(:,2))+max_range/2,...
    mean(Z_real(:,3))-max_range/2,mean(Z_real(:,3))+max_range/2]);
title(title_text);

end