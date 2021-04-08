% This script creates a list of interest points in the mesh. If mesh.matchFixed is true, the mesh will be slightly transformed to include these points

mesh.x.fixPoints = 0;
mesh.y.fixPoints = 0;
mesh.z.fixPoints = [];

if isfield(flowType,'cav')
    for i = 1:length(flowType.cav)
        mesh.x.fixPoints = [mesh.x.fixPoints flowType.cav{i}.x];
        mesh.y.fixPoints = [mesh.y.fixPoints flowType.cav{i}.y];
        mesh.z.fixPoints = [mesh.z.fixPoints flowType.cav{i}.z];
    end
end
if isfield(flowType,'rug')
    for i = 1:length(flowType.rug)
        mesh.x.fixPoints = [mesh.x.fixPoints flowType.rug{i}.x];
        mesh.y.fixPoints = [mesh.y.fixPoints flowType.rug{i}.y];
        mesh.z.fixPoints = [mesh.z.fixPoints flowType.rug{i}.z];
    end
end
if isfield(flowType,'disturb')
    for i = 1:length(flowType.disturb)
        if isfield(flowType.disturb{i},'fitPoints') && flowType.disturb{i}.fitPoints
            mesh.x.fixPoints = [mesh.x.fixPoints flowType.disturb{i}.x];
            mesh.y.fixPoints = [mesh.y.fixPoints flowType.disturb{i}.y];
            mesh.z.fixPoints = [mesh.z.fixPoints flowType.disturb{i}.z];
        else
            if isfield(flowType.disturb{i},'fitPointsX') && flowType.disturb{i}.fitPointsX
                mesh.x.fixPoints = [mesh.x.fixPoints flowType.disturb{i}.x];
            end
            if isfield(flowType.disturb{i},'fitPointsY') && flowType.disturb{i}.fitPointsY
                mesh.y.fixPoints = [mesh.y.fixPoints flowType.disturb{i}.y];
            end
            if isfield(flowType.disturb{i},'fitPointsZ') && flowType.disturb{i}.fitPointsZ
                mesh.z.fixPoints = [mesh.z.fixPoints flowType.disturb{i}.z];
            end
        end
    end
end
if isfield(mesh,'trackedPoints') && isfield(mesh,'fitTrackedPoints') && mesh.fitTrackedPoints
    for i = 1:size(mesh.trackedPoints,1)
        mesh.x.fixPoints = [mesh.x.fixPoints mesh.trackedPoints(i,1)];
        mesh.y.fixPoints = [mesh.y.fixPoints mesh.trackedPoints(i,2)];
        mesh.z.fixPoints = [mesh.z.fixPoints mesh.trackedPoints(i,3)];
    end
end
clear i

mesh.x.fixPoints = unique(mesh.x.fixPoints);
mesh.y.fixPoints = unique(mesh.y.fixPoints);
mesh.z.fixPoints = unique(mesh.z.fixPoints);

mesh.x.fixPoints(isinf(mesh.x.fixPoints) | mesh.x.fixPoints < domain.xi | mesh.x.fixPoints > domain.xf) = [];
mesh.y.fixPoints(isinf(mesh.y.fixPoints) | mesh.y.fixPoints < domain.yi | mesh.y.fixPoints > domain.yf) = [];
mesh.z.fixPoints(isinf(mesh.z.fixPoints) | mesh.z.fixPoints < domain.zi | mesh.z.fixPoints > domain.zf) = [];
