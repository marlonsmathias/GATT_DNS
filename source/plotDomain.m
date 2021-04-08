% This script plots the domain before runtime

colors = lines(6);

figure
hold on

for i = 1:6
    switch i
        case 1
            wallToPlot = boundary.wall.front;
        case 2
            wallToPlot = boundary.wall.back;
        case 3
            wallToPlot = boundary.wall.up;
        case 4
            wallToPlot = boundary.wall.down;
        case 5
            wallToPlot = boundary.wall.right;
        case 6
            wallToPlot = boundary.wall.left;
    end
    
    if ~isempty(wallToPlot)
        switch i
            case {1,2}
                X = mesh.X(wallToPlot(:,[1 1 1 1]));
                Y = mesh.Y(wallToPlot(:,[3 4 4 3]));
                Z = mesh.Z(wallToPlot(:,[5 5 6 6]));
            case {3,4}
                X = mesh.X(wallToPlot(:,[1 2 2 1]));
                Y = mesh.Y(wallToPlot(:,[3 3 3 3]));
                Z = mesh.Z(wallToPlot(:,[5 5 6 6]));
            case {5,6}
                X = mesh.X(wallToPlot(:,[1 2 2 1]));
                Y = mesh.Y(wallToPlot(:,[3 3 4 4]));
                Z = mesh.Z(wallToPlot(:,[5 5 5 5]));
        end

        for j = 1:size(X,1)
            patch(X(j,:),Z(j,:),Y(j,:),colors(i,:))
        end
    end
    
end

corners = boundary.corners.limits;
cornerDir = boundary.corners.dir;
for i = 1:size(corners,1)
    
    if sum(abs(cornerDir(i,:))) == 2
        plot3(mesh.X(corners(i,[1 2])), mesh.Z(corners(i,[5 6])), mesh.Y(corners(i,[3 4])),'r','lineWidth',2);
    else
        plot3(mesh.X(corners(i,1)), mesh.Z(corners(i,5)), mesh.Y(corners(i,3)),'ro');
    end
    
    X = mean(mesh.X(corners(i,[1 2])));
    Y = mean(mesh.Y(corners(i,[3 4])));
    Z = mean(mesh.Z(corners(i,[5 6])));
    
    scale = 0.1;
    
    plot3([X X+cornerDir(i,1)*scale], [Z Z+cornerDir(i,3)*scale], [Y Y+cornerDir(i,2)*scale],'g','lineWidth',2);
    
end

axis tight
axis equal
view(30,30)
