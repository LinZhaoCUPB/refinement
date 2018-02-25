function G = refine_wellregion(G, refine)
% Refine well region by using the radial grids
%
% PARAMETERS:
% G         - global grid structure
% refine    - struct with dimension of (num_of_wellregion * 1) 
%           - refine(well).center - well region/wellbore center coordinates
%           - refine(well).box    - cell index of the refine region:
%                                   [min(x_ind), max(x_ind);
%                                    min(y_ind), max(y_ind)];
%           - refine(well).rW     -  wellbore radius
%           - refine(well).nr     -  number of cells in radial direction
%           - refine(well).rm     -  max radius of the radial grid
%
%
% RETURNS:
% G         - grid structure for global region and well regions
%
% SEE ALSO:
%   `grid_structure`, `tensorGrid`, `radialGrid`


nx = G.cartDims(1);
floc = @(i, j) nx*(j-1) + i;
for k = 1:numel(refine)
    
    % cells in refined region (will be removed)
    boxX = refine(k).box(1,1):refine(k).box(1,2);
    boxY = refine(k).box(2,1):refine(k).box(2,2);
    [boxX, boxY] = meshgrid(boxX, boxY);
    refine_cell = floc(boxX(:), boxY(:));
    refine_cell = find( ismember(G.cells.indexMap, refine_cell) );
    refine(k).refine_cell = refine_cell;
    
    % find boundary of the refined region
    % all nodes in the refined region
    coords = arrayfun(@(c)coords_of_cell(G, c), refine_cell, 'UniformOutput', false);
    coords = cell2mat(coords);
    xmax = max(coords(:,1)); xmin = min(coords(:,1));
    ymax = max(coords(:,2)); ymin = min(coords(:,2));
    refine(k).boundary = [xmin, xmax; ymin, ymax];
    
end
refine_cells = vertcat(refine.refine_cell);

% remove cells in global grid
G = removeCells(G, refine_cells);

% boundary faces (after removeCells)
Gc = computeGeometry(G);
face_cen = Gc.faces.centroids;

for k = 1:numel(refine)
    
    % find boundary faces of well region
    xmin = refine(k).boundary(1,1); xmax = refine(k).boundary(1,2);
    ymin = refine(k).boundary(2,1); ymax = refine(k).boundary(2,2);
    bd_f = any(Gc.faces.neighbors ==0, 2) ...
        & face_cen(:,1) >= xmin - eps & face_cen(:,1) <= xmax + eps...
        & face_cen(:,2) >= ymin - eps & face_cen(:,2) <= ymax + eps;
    bd_f = find(bd_f);
    
    %  find all nodes of boundary faces, to define angles of radial grid
    bd_n    = arrayfun(@(f)nodes_of_face(G, f), bd_f, 'UniformOutput', false);
    bd_n    = cell2mat(bd_n);
    bd_n    = unique(bd_n);
    coords  = G.nodes.coords(bd_n,:);
    center  = refine(k).center;
    coords0 = bsxfun(@minus, coords, center); % for computing angles
    
    % angles of each nodes
    al = cellfun(@compute_angle, num2cell(coords0,2), 'UniformOutput', false);
    al = cell2mat(al);
    [al, ii] = sort(al);
    % sort boundary nodes
    bd_n     = bd_n(ii);
    % sort boundary faces
    bd_f     = sort_bd_faces(G, bd_f, bd_n);
    
    % radial grid parameters
    rm  = refine(k).rm;
    nr  = refine(k).nr;
    rW  = refine(k).rW;
    
    % radial coordinates of nodes
    rr = exp( (linspace(log(rW), log(rm), nr))' );
    
    % define a radial grid
    g = radialGrid(rr, al, center);
    
    % cell number before add
    cell_num0 = G.cells.num;
    % face number before add
    face_num0 = G.faces.num;
    
    % embed the radial grid into the global grid
    G = embed_radialGrid(G, g, bd_f, bd_n);
    
    % cell number after add
    cell_num1 = G.cells.num;
    % face number after add
    face_num1 = G.faces.num;
    
    % cells in refined region
    refine(k).cells = ( cell_num0+1 : cell_num1 )';
    % faces in refined region
    refine(k).faces = ( face_num0+1 : face_num1 )';
    
    na = g.cartDims(2);
    % radial in refined region
    refine(k).radial_cells = ( cell_num0+1 : cell_num1-na )';
    
    % well cells in refined region
    refine(k).well_cells = ( cell_num0+1 : cell_num0+na )';
    
    % wellbore faces in refined region
    intInx = any(G.faces.neighbors(refine(k).faces) == 0, 2);
    refine(k).well_faces = refine(k).faces(intInx);
    
end
G.refine = refine;
G.cells = rmfield(G.cells, 'indexMap');
end

%%
function coords = coords_of_cell(G, c)
% coordinates of cell
ixc = G.cells.facePos;
ixf = G.faces.nodePos;
f = G.cells.faces(ixc(c) : ixc(c + 1) - 1, 1);
e = mcolon(ixf(f), ixf(f + 1) - 1);
nodes  = unique(G.faces.nodes(e, 1));
coords = G.nodes.coords(nodes,:);
end

%%
function nodes = nodes_of_face(G, f)
% nodes of face
ixf    = G.faces.nodePos;
nodes  = G.faces.nodes(ixf(f) : ixf(f + 1) - 1, 1);
end

%%
function al = compute_angle(coord)
get_angle = @(x,y)acos(dot(x,y)/norm(x,2)/norm(y,2));
al = get_angle(coord, [1 0]);
if coord(2) < 0
    al = 2*pi - al;
end
end
%%

function bd_f_sort = sort_bd_faces(G, bd_f, bd_n)
bd_f_sort = zeros(size(bd_f));
% face centroids
pair_nodes = [bd_n(1:end-1), bd_n(2:end)];
pair_nodes = [pair_nodes; [bd_n(end), bd_n(1)]];
for k = 1:length(bd_f)
    nodes = nodes_of_face(G, bd_f(k))';
    ind   = ismember(pair_nodes, nodes);
    ind   =  prod(ind,2)==1 ;
    bd_f_sort(ind) = bd_f(k);
end
end

function G = embed_radialGrid(G, g, bd_f, bd_n)
%% 1. embed radial grid into global grid
cellnum0  = G.cells.num;
facesnum0 = G.faces.num;
nodenum0  = G.nodes.num;

% G.cells -----------------------------------------------------------------
G.cells.num     = cellnum0 + g.cells.num;
G.cells.facePos = [G.cells.facePos; G.cells.facePos(end)-1 + g.cells.facePos(2:end)];
cellfaces = g.cells.faces;
cellfaces(:,1) = cellfaces(:,1) + facesnum0;
G.cells.faces   = [G.cells.faces; cellfaces];

% G.faces -----------------------------------------------------------------
G.faces.num = facesnum0 + g.faces.num;
G.faces.nodePos = [G.faces.nodePos; G.faces.nodePos(end)-1 + g.faces.nodePos(2:end)];

neighbors = g.faces.neighbors;
neighbors(neighbors~=0) = neighbors(neighbors~=0) + cellnum0;
G.faces.neighbors = [G.faces.neighbors; neighbors];

facenodes = g.faces.nodes + nodenum0;
G.faces.nodes = [G.faces.nodes; facenodes];

% G.nodes -----------------------------------------------------------------
G.nodes.num = nodenum0 + g.nodes.num;
G.nodes.coords = [G.nodes.coords; g.nodes.coords];

%% 2. connect the radial grid with global grid
na = g.cartDims(2);
% G.cells -----------------------------------------------------------------
newcells      = G.cells.num + (1:na)';
G.cells.num   = G.cells.num + na;

facePos = cumsum( [1; 4*ones(na,1)] );
G.cells.facePos = [G.cells.facePos; G.cells.facePos(end)-1 + facePos(2:end)];

newfaces  = G.faces.num      + (1:na)';
r0faces   = G.faces.num - na + (1:na)';
r1faces   = bd_f;
faces_a   = [newfaces(1:end-1), newfaces(2:end)];
faces_a   = [faces_a; [newfaces(end), newfaces(1)] ];
cellfaces = [faces_a, r0faces, r1faces];
faceind   = repmat((1:4), na, 1);
cellfaces = [reshape(cellfaces', [], 1), reshape(faceind', [], 1)];

G.cells.faces   = [G.cells.faces; cellfaces];

% G.faces -----------------------------------------------------------------
G.faces.num = G.faces.num + na;

nodePos = cumsum( [1; 2*ones(na,1)] );
G.faces.nodePos = [G.faces.nodePos; G.faces.nodePos(end)-1 + nodePos(2:end)];

neighbors0 = G.faces.neighbors(r0faces,:);
neighbors0(neighbors0 == 0) = newcells;
G.faces.neighbors(r0faces,:) = neighbors0;

neighbors1 = G.faces.neighbors(r1faces,:);
neighbors1(neighbors1 == 0) = newcells;
G.faces.neighbors(r1faces,:) = neighbors1;

neighborsn = [newcells(1:end-1), newcells(2:end)];
neighborsn = [ [newcells(end), newcells(1)]; neighborsn];
G.faces.neighbors = [G.faces.neighbors; neighborsn];

r0nodes   = G.nodes.num - na + (1:na)';
r1nodes   = bd_n;
facenodes = [r0nodes, r1nodes];
facenodes = reshape(facenodes', [], 1);
G.faces.nodes = [G.faces.nodes; facenodes];
end

