function g = radialGrid(rr, al, center)
% Generate a standard grid structure for an individual radial grid.
%
% PARAMETERS:
% rr     - radial  coordinates of nodes, dimension: nr * 1;
% al     - angular coordinates of nodes, dimension: na * 1;
% center - wellbore center coordinates
%
% RETURNS:
% g      - grid structure of mrst
%
% EXAMPLE:
%   % Make a 10 * 30 radial grid
%      nr = 10; rW = 0.2; rm = 1;
%      rr = exp( (linspace(log(rW), log(rm), nr))' );
%      na  = 30;
%      al = linspace(0, 2*pi, na)';
%      al = al(1:end-1);
%      center = [0 0];
%      g = radialGrid(rr, al, center)
%
%   % Plot the grid
%      f = plotGrid(g); axis equal off;
%
% SEE ALSO:
%   `grid_structure`, `tensorGrid`, `cartGrid`


nr = length(rr) - 1;
na = length(al);

% g.cells
% cell numbering: angular (anticlockwise) first, radial second
g.cells.num = nr*na;
g.cells.facePos  = cumsum( [1; 4*ones(g.cells.num,1)] );
g.cells.indexMap  = (1:g.cells.num)';
faces = arrayfun(@(cell)faces_of_cell(cell, nr, na), g.cells.indexMap, 'UniformOutput', false);
g.cells.faces     = cell2mat(faces);

% g.faces
% face numbering: angular (anticlockwise) first, radial second
% (1:na*nr), angular faces; (na*nr)+ (1:na*(nr+1)), angular faces, 
g.faces.num    = (na*nr) + na*(nr+1);
g.faces.nodePos = cumsum( [1; 2*ones(g.faces.num,1)] );
nodes_a = arrayfun(@(face)nodes_of_face_angular(face, nr, na), (1:na*nr)',              'UniformOutput', false);
nodes_r = arrayfun(@(face)nodes_of_face_radial (face, nr, na), (na*nr)+ (1:na*(nr+1))', 'UniformOutput', false);
g.faces.nodes   = cell2mat( vertcat(nodes_a, nodes_r) );

neighbors_a = arrayfun(@(face)neighbors_of_face_angular(face, nr, na), (1:na*nr)',              'UniformOutput', false);
neighbors_r = arrayfun(@(face)neighbors_of_face_radial (face, nr, na), (na*nr)+ (1:na*(nr+1))', 'UniformOutput', false);
g.faces.neighbors   = cell2mat( vertcat(neighbors_a, neighbors_r) );


%g.nodes
g.nodes.num    = na*(nr+1);
coords_x = rr * cos(al') + center(1);
coords_y = rr * sin(al') + center(2);
g.nodes.coords = [reshape(coords_x', [], 1), reshape(coords_y', [], 1)];

% others
g.cartDims = [nr, na];
g.type     = {'radialGrid'};
g.griddim  = 2;
g.center   = center;
end


%%
function faces = faces_of_cell(cell, nr, na)
% define faces of each cell
% face numbering:angular (anticlockwise) first, radial second
% 1 - angular- direction
% 2 - angular+ direction
% 3 - radial-  direction
% 4 - radial+  direction
faces = zeros(4,2);
faces(1,1) = cell;
faces(2,1) = cell+1;
faces(3,1) = nr*na + cell;
faces(4,1) = nr*na + cell+na;

% joint the first cell and the last cell in each angular loop
if rem( cell/na, 1 ) == 0
    faces(2,1) = cell - na + 1;
end

faces(:,2) = (1:4)';
end

%%
function neighbors = neighbors_of_face_angular(face, nr, na)
% define neighbors of each angular face
% neighbors(1) - angular-
% neighbors(2) - angular+
neighbors = zeros(1,2);
neighbors(1) = face-1;
neighbors(2) = face;

if rem( (face-1)/na, 1 ) == 0
    neighbors(1) = face + na - 1;
end
% in accordance with radial
nr = nr;
end

%%
function neighbors = neighbors_of_face_radial(face, nr, na)
% define neighbors of each radial face
% neighbors(1) - radial-
% neighbors(2) - radial+
neighbors = zeros(1,2);
face0 = face - nr*na;
neighbors(1) = face0 - na;
neighbors(2) = face0;

% cells (1 ~ nr*na)
neighbors( neighbors < 1       ) = 0;
neighbors( neighbors > nr*na   ) = 0;
end

%%
function nodes = nodes_of_face_angular(face, nr, na)
% define nodes of each angular face
% node numbering:angular (anticlockwise) first, radial second
nodes = zeros(2,1);
nodes(1) = face;
nodes(2) = face + na;
% in accordance with radial
nr = nr;
end

%%
function nodes = nodes_of_face_radial(face, nr, na)
% define nodes of each radial face
% node numbering:angular (anticlockwise) first, radial second
nodes = zeros(2,1);
face0 = face - nr*na;
nodes(1) = face0;
nodes(2) = face0 + 1;

if rem( face0/na, 1 ) == 0
    nodes(2) = face0 - na + 1;
end
end