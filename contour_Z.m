function contour_Z(x,y,z,lambda,color,psize)
%%Input x,y,z data which are randomly spaced
%   Plot a contour from those random data
%   First, use meshgrid to generate regularly spaced coordinates, and second 
%       use the griddata to generate interpolation between those data
%       pionts, third plot the contour
%   x: x coordinates
%   y: y coordinates
%   z: height for x,y (could be multiple columns; one plot for each column

xlin = linspace(min(x),max(x),33);
ylin = linspace(min(y),max(y),33);
[X,Y] = meshgrid(xlin,ylin);

np = size(z,2);

figure();
hold on;
for i = 1:np
    Z = griddata(x,y,z(:,i),X,Y,'cubic');
    subplot(np,2,2*i-1);
    contour(X,Y,Z,'ShowText','on');
    title(['Contour plot of ',num2str(i),' coordinates, \lambda=',num2str(lambda(i))]);
    xlabel('Z_1');
    ylabel('Z_2');
    
    subplot(np,2,2*i);
    mesh(X,Y,Z);
    title(['Surface plot of ',num2str(i),' coordinates']);
    hold on;
    axis tight;
    scatter3(x,y,z(:,i),psize,color,'filled')
    xlabel('Z_1');
    ylabel('Z_2');
    zlabel(['Z_{est}']);
    %view(az,el);
end


end