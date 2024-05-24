%% Subroutine for a periodic box
% To be called from getBoundaryConditions

flowRegion = true(mesh.nx,mesh.ny,mesh.nz);

findWallsForBoundaries
