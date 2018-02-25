%% Well Region Refinement by Radial Grid
clear
mrstModule add ad-core ad-blackoil ad-props mrst-gui
close all
%% Reservoir geometry
% only 2D grid is supported
G = cartGrid([30 30], [600, 600]*meter);
%% Define well regions
dxy = [20, 20]; % 600/30
refine = struct('center',[], 'box', []);

% well region 1 -----------------------------------------------------------
% centeral cell index
cc_ind1 = [5, 10];
% well region/wellbore center coordinates, center = [x_coord, y_coord];
refine(1).center = dxy.*(cc_ind1 - 1) + dxy/2;
% cell index of the refine region
% [min(x_ind), max(x_ind); min(y_ind), max(y_ind)];
sx = 3; sy = 3;
refine(1).box    = [cc_ind1(1) - sx, cc_ind1(1) + sx; cc_ind1(2) - sy, cc_ind1(2) + sy];
% wellbore radius
refine(1).rW     = 0.2;
% number of cells in radial direction
refine(1).nr     = 10;
% max radius of the radial grid
refine(1).rm     = 50;

% well region 2 -----------------------------------------------------------
cc_ind2 = [25, 5];
refine(2).center = dxy.*(cc_ind2 - 1) + dxy/2;
sx = 3; sy = 3;
refine(2).box    = [cc_ind2(1) - sx, cc_ind2(1) + sx; cc_ind2(2) - sy, cc_ind2(2) + sy];
refine(2).rW     = 0.2;
refine(2).nr     = 10;
refine(2).rm     = 50;

% well region 3 -----------------------------------------------------------
cc_ind3 = [25, 25];
refine(3).center = dxy.*(cc_ind3 - 1) + dxy/2;
sx = 3; sy = 2;
refine(3).box    = [cc_ind3(1) - sx, cc_ind3(1) + sx; cc_ind3(2) - sy, cc_ind3(2) + sy];
refine(3).rW     = 0.2;
refine(3).nr     = 10;
refine(3).rm     = 40;

nx = G.cartDims(1);
floc = @(i, j) nx*(j-1) + i;

figure; hold on
plotGrid(G)
col = {'r','g','b'};
for k = 1:numel(refine)
    boxX = refine(k).box(1,1):refine(k).box(1,2);
    boxY = refine(k).box(2,1):refine(k).box(2,2);
    [boxX, boxY] = meshgrid(boxX, boxY);
    refine_cell = floc(boxX(:), boxY(:));
    refine_cell = find( ismember(G.cells.indexMap, refine_cell) );
    plotGrid(G, refine_cell, 'facecol', col{k})
end
axis equal off
legend({'global','refinement1','refinement2','refinement3'}, 'Location', 'NorthWest')
%% Refine well regions
% 1. grid structure for refined grid
G_R = refine_wellregion(G, refine);
G_R = computeGeometry(G_R);

figure;
subplot(2,2,1);plotGrid(G_R);axis equal off; title('global');
for k = 1:3
    subplot(2,2,k+1);plotGrid(G_R, G_R.refine(k).cells);axis equal off 
    title(sprintf('refinement %d', k));
end
drawnow

% 2. grid structure for cart grid
G   = computeGeometry(G);
%% Set up rock properties
% 1. rock for refined grid
nc_R = G_R.cells.num;
rock_R.perm = 500*milli*darcy * ones(nc_R,2);
rock_R.poro = 0.3             * ones(nc_R,1);

% 2. rock for cart grid
nc   = G.cells.num;
rock.perm   = 500*milli*darcy * ones(nc,2);
rock.poro   = 0.3             * ones(nc,1);
%% Set up simulation model
fluid = initSimpleADIFluid('mu',    [1, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
                       
% Constant oil compressibility
c        = 0.001/barsa;
p_ref    = 300*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

gravity reset on

% 1. model for refined grid
model_R = TwoPhaseOilWaterModel(G_R, rock_R, fluid);
% modify transmissibility and pv based on radial formulae
model_R.operators = operators_radialGrid(G_R, rock_R, model_R.operators);

% 2. model for cart grid
model = TwoPhaseOilWaterModel(G, rock, fluid);
%% Define wells and simulation schedule
% 1. refined grid, the wellbore is shown explicitly, the well cells are the
% surrounding cells, and the well index equals wellbore face transmissibility
% add injector
W_R = [];
for k = 1
    wc = G_R.refine(k).well_cells;
    wf = G_R.refine(k).well_faces;
    % well index = wellbore face transmissibility
    WI = model_R.operators.T_all(wf);
    W_R = addWell(W_R, G_R, rock_R, wc, 'Name', 'I1', ...
     'Val', 400*barsa, 'Type', 'bhp', 'sign', 1, 'Comp_i', [1 0], 'WI', WI);
end
% add producer
for k = 2:3
    wc = G_R.refine(k).well_cells;
    wf = G_R.refine(k).well_faces;
    % well index = wellbore face transmissibility
    WI = model_R.operators.T_all(wf);
    W_R = addWell(W_R, G_R, rock_R, wc, 'Name', sprintf('P%d', k-1), ...
    'Val', 150*barsa, 'Type', 'bhp', 'sign', -1, 'Comp_i', [1 1], 'WI', WI);
end

% 2. cart grid, the well is virtual, and the well index is calculated by
% PI = 2*pi*k*h/ln(req/rw), req = 0.14*sqrt(dx^2 + dy^2) (Peaceman, 1983)
% add injector
W = [];
W = verticalWell(W, G, rock,  cc_ind1(1), cc_ind1(2), 1, 'Name', 'I1', ...
    'Val', 400*barsa, 'Type', 'bhp', 'sign',  1, 'Comp_i', [1 0]);

% add producer
W = verticalWell(W, G, rock,  cc_ind2(1), cc_ind2(2), 1, 'Name', 'P1', ...
    'Val', 150*barsa, 'Type', 'bhp', 'sign', -1, 'Comp_i', [1 1]);

W = verticalWell(W, G, rock,  cc_ind3(1), cc_ind3(2), 1, 'Name', 'P2', ...
    'Val', 150*barsa, 'Type', 'bhp', 'sign', -1, 'Comp_i', [1 1]);


% Set up time steps. Since the grid near the wellbore is very fine, we define
% small time steps in the beginning of simulation
timesteps1 = exp( (linspace(log(0.1*hour), log(10*day), 20))' );
timesteps2 = ones(100,1)*10*day;
timesteps = [timesteps1; timesteps2];

% Set up the schedule containing both the wells and the timesteps
schedule_R = simpleSchedule(timesteps, 'W', W_R);
schedule   = simpleSchedule(timesteps, 'W', W  );
%% define initial state
state0_R = initResSol(G_R, 300*barsa, [0 1]);
state0   = initResSol(G,   300*barsa, [0 1]);
%% Simulate base case
[wellSols_R, states_R, report_R] = simulateScheduleAD(state0_R, model_R, schedule_R);

[wellSols, states, report] = simulateScheduleAD(state0, model, schedule);
%% Compare the simluation results
% water rate
qWs_R = cellfun(@(x)horzcat(x.qWs) * day, wellSols_R , 'UniformOutput', false);
qWs   = cellfun(@(x)horzcat(x.qWs) * day, wellSols,    'UniformOutput', false);
qWs_R = cell2mat(qWs_R);  qWs   = cell2mat(qWs);

% oil rate
qOs_R = cellfun(@(x)horzcat(x.qOs) * day, wellSols_R , 'UniformOutput', false);
qOs   = cellfun(@(x)horzcat(x.qOs) * day, wellSols,    'UniformOutput', false);
qOs_R = cell2mat(qOs_R);  qOs   = cell2mat(qOs);

RT_R = report_R.ReservoirTime/day;
RT   = report.ReservoirTime  /day;

figure;hold on; grid on; box on
plot(RT_R, qWs_R(:,1), 'r.-')
plot(RT,   qWs  (:,1), 'g.-')
xlim([0 1200])
xlabel('Simulation Time (day)')
ylabel('Injetion Rate (m^3/day)')
title('Injetion Rate of I1')
legend({'refined grid', 'cartGrid'}, 'Location','NorthWest')


figure;hold on; grid on; box on
plot(RT_R, -qWs_R(:,2), 'r.-')
plot(RT,   -qWs  (:,2), 'g.-')
plot(RT_R, -qOs_R(:,2), 'b.-')
plot(RT,   -qOs  (:,2), 'm.-')
xlim([0 1200])
xlabel('Simulation Time (day)')
ylabel('Production Rate (m^3/day)')
title('Production Rate of P1')
legend({'Water Rate - refined grid', 'Water Rate - cartGrid',...
    'Oil Rate - refined grid', 'Oil Rate - cartGrid'}, 'Location','NorthWest')

figure;hold on; grid on; box on
plot(RT_R, -qWs_R(:,3), 'r.-')
plot(RT,   -qWs  (:,3), 'g.-')
plot(RT_R, -qOs_R(:,3), 'b.-')
plot(RT,   -qOs  (:,3), 'm.-')
xlim([0 1200])
xlabel('Simulation Time (day)')
ylabel('Production Rate (m^3/day)')
title('Production Rate of P2')
legend({'Water Rate - refined grid', 'Water Rate - cartGrid',...
    'Oil Rate - refined grid', 'Oil Rate - cartGrid'}, 'Location','NorthWest')


figure;
subplot(2,2,1); axis equal off; title('oil saturation, refined grid')
plotCellData(G_R, states_R{end}.s(:,2));colorbar

subplot(2,2,2); axis equal off; title('oil saturation, cart grid')
plotCellData(G,   states{end}.s(:,2));colorbar

subplot(2,2,3); axis equal off; title('pressure (barsa), refined grid')
plotCellData(G_R, states_R{end}.pressure/barsa);colorbar

subplot(2,2,4); axis equal off; title('pressure (barsa), cart grid')
plotCellData(G,   states{end}.pressure/barsa);colorbar
