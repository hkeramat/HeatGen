% File: create_params.m
clear all;
% Import COMSOL classes
import com.comsol.model.*
import com.comsol.model.util.*
% (2) Create or load a COMSOL model
modelTags = ModelUtil.tags();
% Remove each model using its tag
for i = 1:length(modelTags)
    ModelUtil.remove(modelTags(i));
end

% Confirm all models are removed
%disp(['Number of models remaining: ', num2str(ModelUtil.size())]);
model = ModelUtil.create('Model1');%mphload('/home/morteza/Desktop/two_shape_merged_bezier.mph');

% Define parameters
model.param.set('H', '1[m]');
model.param.set('cp_f', '1006[J/(kg*K)]');
model.param.set('mu_f', '1.94E-5[Pa*s]');
model.param.set('rho_f', '1.204[kg/m^3]');
model.param.set('k_f', '0.024[W/(m*K)]');
model.param.set('rho_s', '2700[kg/m^3]');
model.param.set('k_s', '237[W/(m*K)]');
model.param.set('cp_s', '900[J/(kg*K)]');
model.param.set('T_inlet', '293.15[K]');
model.param.set('H_w_dim', '1[m]');
model.param.set('u_inlet', '1[m/s]');
model.param.set('Q0', '1E5[W/m^2]');
model.param.set('fluid_area_width', '0.25*H_dim');
model.param.set('fluid_area_height', 'H_dim');
model.param.set('fluid_solid_area_height', 'H_dim');
model.param.set('fluid_solid_area_width', 'H_dim');
model.param.set('solid_area_width', 'fluid_solid_area_width');
model.param.set('solid_area_height', 'fluid_solid_area_height');
model.param.set('x_center_fluid_area', '(fluid_area_width+fluid_solid_area_width)/2');
model.param.set('y_center_fluid_area', 'fluid_area_height/2');
model.param.set('x_center_fluid2_area', '-(fluid_area_width+fluid_solid_area_width)/2');
model.param.set('y_center_fluid2_area', 'fluid_area_height/2');
model.param.set('x_center_fluid_solid_area', '0[m]/x0');
model.param.set('y_center_fluid_solid_area', 'fluid_solid_area_height/2');
model.param.set('x0', '10[mm]');
model.param.set('H_dim', '1[m]');
model.param.set('y_displacement', '1.5*H_dim');
model.param.set('x_center_solid_area', 'x_center_fluid_solid_area');
model.param.set('y_center_solid_area', 'y_center_fluid_solid_area - y_displacement');
model.param.set('dz_bp', 'H_dim/8');
model.param.set('dz_channel', '1.5*H_dim');
model.param.set('h_f', '50[W/(m^2*K)]');
model.param.set('h_s', '48000[W/(m^2*K)]');
% Additional parameters
% param_list = {'x21', '0.24*2-0.2'; 'x22', '0.2*2-0.2'; 'x23', '0.14*2-0.2'; 'x24', '0.08*2-0.2';
%               'x25', '0.04*2-0.2'; 'x26', '0.03*2-0.2'; 'x27', '0.03*2-0.2'; 'x28', '0.05*2-0.2';
%               'x29', '0.08*2-0.2'; 'x210', '0.12*2-0.2'; 'x211', '0.16*2-0.2'; 'x212', '0.2*2-0.2';
%               'y21', '0+0.2+0.5'; 'y22', '0.02*2+0.2+0.5'; 'y23', '0.02*2+0.2+0.5'; 'y24', '0.02*2+0.2+0.5';
%               'y25', '0.02*2+0.2+0.5'; 'y26', '0.01*2+0.2+0.5'; 'y27', '0+0.2+0.5'; 'y28', '-0.01*2+0.2+0.5';
%               'y29', '-0.01*2+0.2+0.5'; 'y210', '-0.01*2+0.2+0.5'; 'y211', '-0.01*2+0.2+0.5'; 'y212', '-0.01*2+0.2+0.5';
%               'x11', '0.24*2-0.2'; 'x12', '0.2*2-0.2'; 'x13', '0.14*2-0.2'; 'x14', '0.08*2-0.2';
%               'x15', '0.04*2-0.2'; 'x16', '0.03*2-0.2'; 'x17', '0.03*2-0.2'; 'x18', '0.05*2-0.2';
%               'x19', '0.08*2-0.2'; 'x110', '0.12*2-0.2'; 'x111', '0.16*2-0.2'; 'x112', '0.2*2-0.2';
%               'y11', '0+0.2+0.1'; 'y12', '0.02*2+0.2+0.1'; 'y13', '0.02*2+0.2+0.1'; 'y14', '0.02*2+0.2+0.1';
%               'y15', '0.02*2+0.2+0.1'; 'y16', '0.01*2+0.2+0.1'; 'y17', '0+0.2+0.1'; 'y18', '-0.01*2+0.2+0.1';
%               'y19', '-0.01*2+0.2+0.1'; 'y110', '-0.01*2+0.2+0.1'; 'y111', '-0.01*2+0.2+0.1'; 'y112', '-0.01*2+0.2+0.1'};
% 
% for i = 1:size(param_list,1)
%     model.param.set(param_list{i,1}, param_list{i,2});
% end
data = load('pmcopy3_parameters.mat');
fid = fopen('pmcopy3_iteration_count.txt', 'r');
iteration = fscanf(fid, '%d');  % or '%f' if float
fclose(fid);
% Iterate over the fields and set parameters in COMSOL
paramNames = fieldnames(data);
for i = 1:numel(paramNames)
    paramName = paramNames{i};
    paramValue = data.(paramName);
    model.param.set(paramName, num2str(paramValue), '');
end



x0=mphevaluate(model,'x0');
fluid_solid_area_width=mphevaluate(model,'fluid_solid_area_width');
fluid_solid_area_height=mphevaluate(model,'fluid_solid_area_height');
x_center_fluid_solid_area=mphevaluate(model,'x_center_fluid_solid_area');
y_center_fluid_solid_area=mphevaluate(model,'y_center_fluid_solid_area');
T_inlet=mphevaluate(model,'T_inlet');
rho_f=mphevaluate(model,'rho_f');
H=mphevaluate(model,'H');
dz_channel=mphevaluate(model,'dz_channel');
u_inlet=mphevaluate(model,'u_inlet');
Q0=mphevaluate(model,'Q0');
cp_f=mphevaluate(model,'cp_f');
num_fins=2;





model.component().create("comp1");
model.geom().create("geom1", 2);



% Rectangle r1
model.geom('geom1').create('r1', 'Rectangle');
model.geom('geom1').feature('r1').set('base', 'center');
model.geom('geom1').feature('r1').set('pos', {'x_center_fluid_area' 'y_center_fluid_area'});
model.geom('geom1').feature('r1').set('size', {'fluid_area_width' 'fluid_area_height'});
model.geom('geom1').feature('r1').set('createselection', 'on');
model.geom('geom1').feature('r1').label('Rectangle R1');

% Rectangle r2
model.geom('geom1').create('r2', 'Rectangle');
model.geom('geom1').feature('r2').set('base', 'center');
model.geom('geom1').feature('r2').set('pos', {'x_center_fluid2_area' 'y_center_fluid2_area'});
model.geom('geom1').feature('r2').set('size', {'fluid_area_width' 'fluid_area_height'});
model.geom('geom1').feature('r2').set('createselection', 'on');
model.geom('geom1').feature('r2').label('Rectangle R2');

% Rectangle r3
model.geom('geom1').create('r3', 'Rectangle');
model.geom('geom1').feature('r3').set('base', 'center');
model.geom('geom1').feature('r3').set('pos', {'x_center_fluid_solid_area' 'y_center_fluid_solid_area'});
model.geom('geom1').feature('r3').set('size', {'fluid_solid_area_width' 'fluid_solid_area_height'});
model.geom('geom1').feature('r3').set('createselection', 'on');
model.geom('geom1').feature('r3').label('Rectangle R3');

% Rectangle r4
model.geom('geom1').create('r4', 'Rectangle');
model.geom('geom1').feature('r4').set('base', 'center');
model.geom('geom1').feature('r4').set('pos', {'x_center_solid_area' 'y_center_solid_area'});
model.geom('geom1').feature('r4').set('size', {'solid_area_width' 'solid_area_height'});
model.geom('geom1').feature('r4').set('createselection', 'on');
model.geom('geom1').feature('r4').label('Rectangle R4');


model.component("comp1").geom("geom1").create("cc1", "CompositeCurve");
model.component("comp1").geom("geom1").feature("cc1").create("cb1", "CubicBezier");
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "x11", 0, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "y11", 1, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "x12", 0, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "y12", 1, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "x13", 0, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "y13", 1, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "x14", 0, 3);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "y14", 1, 3);


model.component("comp1").geom("geom1").feature("cc1").create("cb2", "CubicBezier");
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "x14", 0, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "y14", 1, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "x15", 0, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "y15", 1, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "x16", 0, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "y16", 1, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "x17", 0, 3);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "y17", 1, 3);

model.component("comp1").geom("geom1").feature("cc1").create("cb3", "CubicBezier");
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "x17", 0, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "y17", 1, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "y18", 1, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "x18", 0, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "x19", 0, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "y19", 1, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "x110", 0, 3);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "y110", 1, 3);

model.component("comp1").geom("geom1").feature("cc1").create("cb4", "CubicBezier");
model.component("comp1").geom("geom1").feature("cc1").feature("cb4").setIndex("p", "x110", 0, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb4").setIndex("p", "y110", 1, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb4").setIndex("p", "x11", 0, 3);
model.component("comp1").geom("geom1").feature("cc1").feature("cb4").setIndex("p", "y11", 1, 3);
model.component("comp1").geom("geom1").feature("cc1").feature("cb4").setIndex("p", "x111", 0, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb4").setIndex("p", "y111", 1, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb4").setIndex("p", "x112", 0, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb4").setIndex("p", "y112", 1, 2);

model.component("comp1").geom("geom1").feature("cc1").set("createselection", "on");
model.component("comp1").geom("geom1").feature("cc1").set("selresultshow", "all");
model.component("comp1").geom("geom1").feature("cc1").set("type", "solid");

model.component("comp1").geom("geom1").create("cc2", "CompositeCurve");
model.component("comp1").geom("geom1").feature("cc2").create("cb1", "CubicBezier");
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "x21", 0, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "y21", 1, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "x22", 0, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "y22", 1, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "x23", 0, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "y23", 1, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "x24", 0, 3);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "y24", 1, 3);

model.component("comp1").geom("geom1").feature("cc2").create("cb2", "CubicBezier");
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "x24", 0, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "y24", 1, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "x25", 0, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "y25", 1, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "x26", 0, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "y26", 1, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "x27", 0, 3);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "y27", 1, 3);

model.component("comp1").geom("geom1").feature("cc2").create("cb3", "CubicBezier");
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "x27", 0, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "y27", 1, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "y28", 1, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "x28", 0, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "x29", 0, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "y29", 1, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "x210", 0, 3);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "y210", 1, 3);

model.component("comp1").geom("geom1").feature("cc2").create("cb4", "CubicBezier");
model.component("comp1").geom("geom1").feature("cc2").feature("cb4").setIndex("p", "x210", 0, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb4").setIndex("p", "y210", 1, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb4").setIndex("p", "x21", 0, 3);
model.component("comp1").geom("geom1").feature("cc2").feature("cb4").setIndex("p", "y21", 1, 3);
model.component("comp1").geom("geom1").feature("cc2").feature("cb4").setIndex("p", "x211", 0, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb4").setIndex("p", "y211", 1, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb4").setIndex("p", "x212", 0, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb4").setIndex("p", "y212", 1, 2);

model.component("comp1").geom("geom1").feature("cc2").set("createselection", "on");
model.component("comp1").geom("geom1").feature("cc2").set("selresultshow", "all");
model.component("comp1").geom("geom1").feature("cc2").set("type", "solid");








y_displacement=mphevaluate(model,'y_displacement');


geomName = 'geom1';
compName='comp1';
geom = model.component(compName).geom(geomName);




copy1_Tag='copy1';
geom.feature.create(copy1_Tag, 'Copy');
geom.feature(copy1_Tag).selection('input').set({'cc1'}); % Adjust tags
geom.feature(copy1_Tag).set('displ', {'0', '-y_displacement'});
geom.feature(copy1_Tag).set('selresult', true);
model.component("comp1").geom("geom1").feature("copy1").set("propagatesel", false);

copy2_Tag='copy2';
geom.feature.create(copy2_Tag, 'Copy');
geom.feature(copy2_Tag).selection('input').set({'cc2'}); % Adjust tags
geom.feature(copy2_Tag).set('displ', {'0', '-y_displacement'});
geom.feature(copy2_Tag).set('selresult', true);
model.component("comp1").geom("geom1").feature("copy2").set("propagatesel", false);


copy3_Tag='copy3';
geom.feature.create(copy3_Tag, 'Copy');
geom.feature(copy3_Tag).selection('input').set({'cc1'}); % Adjust tags
geom.feature(copy3_Tag).set('displ', {'0', 'y_displacement'});
geom.feature(copy3_Tag).set('selresult', true);
model.component("comp1").geom("geom1").feature("copy3").set("propagatesel", false);


copy4_Tag='copy4';
geom.feature.create(copy4_Tag, 'Copy');
geom.feature(copy4_Tag).selection('input').set({'cc2'}); % Adjust tags
geom.feature(copy4_Tag).set('displ', {'0', 'y_displacement'});
geom.feature(copy4_Tag).set('selresult', true);
model.component("comp1").geom("geom1").feature("copy4").set("propagatesel", false);







fin_in_liquid_Tag='fin_in_liquid';
geom.feature.create(fin_in_liquid_Tag, 'Union');
geom.feature(fin_in_liquid_Tag).selection('input').set({'cc1','cc2'}); % Adjust tags %,'cc3'
feature = geom.feature(fin_in_liquid_Tag);
% Activate Result object selection
feature.set('selresult', true);
% Deactivate Keep interior boundaries
feature.set('intbnd', false);




% Create a difference operation
differenceTag = 'liquid_solid'; % Unique tag for the difference operation
geom.feature.create(differenceTag, 'Difference');
% Specify the object to subtract from (the base object)
geom.feature(differenceTag).selection('input').set({'r3'});
% Specify the objects to subtract
geom.feature(differenceTag).selection('input2').set({'fin_in_liquid'}); 
% Delete the base object (input)
geom.feature(differenceTag).set('keep', true);
% Keep interior boundaries
geom.feature(differenceTag).set('intbnd', true);
feature = geom.feature(differenceTag);
% Activate Result object selection
feature.set('selresult', true);


liquid_tot_Tag='liquid_tot';
geom.feature.create(liquid_tot_Tag, 'Union');
geom.feature(liquid_tot_Tag).selection('input').set({'r1','r2','liquid_solid'}); % Adjust tags
feature = geom.feature(liquid_tot_Tag);
% Activate Result object selection
feature.set('selresult', true);
% Deactivate Keep interior boundaries
feature.set('intbnd', true);
model.component("comp1").geom("geom1").feature("fin_in_liquid").set("selresultshow", "all");


fin_in_solid_Tag='fin_in_solid';
geom.feature.create(fin_in_solid_Tag, 'Union');
geom.feature(fin_in_solid_Tag).selection('input').set({'copy1','copy2'}); % Adjust tags %,'copy3'
feature = geom.feature(fin_in_solid_Tag);
% Activate Result object selection
feature.set('selresult', true);
% Deactivate Keep interior boundaries
feature.set('intbnd', false);


% Use a Difference operation: object1 minus the intersection.
differenceTag = 'solid_liquid';
geom.feature.create(differenceTag, 'Difference');
% The base (input) is 'object1', we subtract the intersection ('int1')
geom.feature(differenceTag).selection('input').set({'r4'});
geom.feature(differenceTag).selection('input2').set({'fin_in_solid'});
feature = geom.feature(differenceTag);
% Activate Result object selection
feature.set('selresult', true);
feature.set('keep', true);
% Deactivate Keep interior boundaries
feature.set('intbnd', true);

feature = geom.feature('r1');
feature.set('selresult', true);

feature = geom.feature('r4');
feature.set('selresult', true);

feature = geom.feature('r2');
feature.set('selresult', true);

geom.run();

domainIndices_fin_in_solid = model.component('comp1').selection('geom1_fin_in_solid_dom').entities();
domainIndices_solid_liquid=model.component('comp1').selection('geom1_solid_liquid_dom').entities();
domainIndices_right_rectangle = model.component('comp1').selection('geom1_r1_dom').entities();
domainIndices_left_rectangle = model.component('comp1').selection('geom1_r2_dom').entities();
domainIndices_mid_rectangle = model.component('comp1').selection('geom1_liquid_solid_dom').entities();
domainIndices_down_rectangle = model.component('comp1').selection('geom1_r4_dom').entities();
domainIndices_fin_in_liquid = model.component('comp1').selection('geom1_fin_in_liquid_dom').entities();
domainIndices_liquid_solid = model.component('comp1').selection('geom1_liquid_solid_dom').entities();
domainIndices_liquid_tot = model.component('comp1').selection('geom1_liquid_tot_dom').entities();
domainIndices_copy3 = model.component('comp1').selection('geom1_copy3_dom').entities();
domainIndices_copy4 = model.component('comp1').selection('geom1_copy4_dom').entities();
bondIndices_fin_in_liquid = model.component('comp1').selection('geom1_fin_in_liquid_bnd').entities();


% Convert all domain index arrays into row vectors
finInSolid = domainIndices_fin_in_solid(:).';          % Make sure it's a row vector
solidLiquid = domainIndices_solid_liquid(:).';         % Row vector
finInLiquid=domainIndices_fin_in_liquid(:).';
liquidSolid=domainIndices_liquid_solid(:).';

allDomains = [finInSolid, solidLiquid];
allDomains2 = [finInLiquid,liquidSolid];


% Get the adjacency matrix for points (vertices) connected to domains
adj = mphgetadj(model, 'geom1', 'point', 'domain');
pointsRect3 = adj{domainIndices_mid_rectangle};
% disp('Points (vertices) of Rectangle 3:');
% disp(pointsRect3);

coordsRect3 = mphgetcoords(model, 'geom1', 'point', pointsRect3);
% disp('Coordinates of points in Rectangle 3:');
% disp(coordsRect3);

% Find the smallest x-coordinate
minX = min(coordsRect3(1, :));
maxX = max(coordsRect3(1, :));
% Get indices of vertices with the lowest x-coordinate
candidateIndices = find(coordsRect3(1, :) == minX);
candidateIndices2 = find(coordsRect3(1, :) == maxX);
% Among these candidates, find the one with the largest y-coordinate
[~, idx] = max(coordsRect3(2, candidateIndices)); % Largest y-coordinate
[~, idx2] = min(coordsRect3(2, candidateIndices)); % Smallest y-coordinate
[~, idx3] = min(coordsRect3(2, candidateIndices2)); % Largest y-coordinate
upperLeftPoint3 = candidateIndices(idx);
downLeftPoint3 = candidateIndices(idx2);
upperrightPoint3 = candidateIndices2(idx3);
% Get the vertex ID for the upper-right corner
upperLeftVertex3 = pointsRect3(upperLeftPoint3);
downLeftVertex3 = pointsRect3(downLeftPoint3);
upperrightVertex3 = pointsRect3(upperrightPoint3);
% disp(['Upper-left corner vertex ID for the rectangle 3: ', num2str(upperLeftVertex3)]);

pointsRect4 = adj{domainIndices_down_rectangle};
% disp('Points (vertices) of Rectangle 4:');
% disp(pointsRect4);

coordsRect4 = mphgetcoords(model, 'geom1', 'point', pointsRect4);
% disp('Coordinates of points in Rectangle 4:');
% disp(coordsRect4);

minX = min(coordsRect4(1, :));
maxX = max(coordsRect4(1, :));
% Get indices of vertices with the lowest x-coordinate
candidateIndices = find(coordsRect4(1, :) == minX);
candidateIndices2 = find(coordsRect4(1, :) == maxX);
% Among these candidates, find the one with the largest y-coordinate
[~, idx] = max(coordsRect4(2, candidateIndices)); % Largest y-coordinate
[~, idx2] = min(coordsRect4(2, candidateIndices)); % Smallest y-coordinate
[~, idx3] = min(coordsRect4(2, candidateIndices2)); % Largest y-coordinate
upperLeftPoint4 = candidateIndices(idx);
downLeftPoint4 = candidateIndices(idx2);
upperrightPoint4 = candidateIndices2(idx3);
% Get the vertex ID for the upper-right corner
upperLeftVertex4 = pointsRect4(upperLeftPoint4);
downLeftVertex4 = pointsRect4(downLeftPoint4);
upperrightVertex4 = pointsRect4(upperrightPoint4);
% disp(['Upper-left corner vertex ID for the rectangle 4: ', num2str(upperLeftVertex4)]);


% Create a named selection for the source domains
 model.component('comp1').selection.create('sourceDomSel', 'Explicit');
 model.component('comp1').selection('sourceDomSel').geom('geom1', 2);  % 2 for 2D domains
 model.component('comp1').selection('sourceDomSel').set(allDomains);  % the identified domain indices

 % Create a named selection for the source domains
 model.component('comp1').selection.create('sourceDomSel2', 'Explicit');
 model.component('comp1').selection('sourceDomSel2').geom('geom1', 2);  % 2 for 2D domains
 model.component('comp1').selection('sourceDomSel2').set(allDomains2);  % the identified domain indices


% Create the LinearExtrusion coupling
model.component('comp1').cpl.create('linext1', 'LinearExtrusion');
model.component('comp1').cpl.create('linext2', 'LinearExtrusion');

 % Assign the domain selection to the coupling
model.component('comp1').cpl('linext1').selection.named('sourceDomSel');
model.component('comp1').cpl('linext2').selection.named('sourceDomSel2');
% Assign the point selections to LinearExtrusion
model.component('comp1').cpl('linext1').selection('srcvertex1').set(upperLeftVertex4);
model.component('comp1').cpl('linext1').selection('dstvertex1').set(upperLeftVertex3);
model.component('comp1').cpl('linext1').selection('srcvertex2').set(downLeftVertex4);
model.component('comp1').cpl('linext1').selection('dstvertex2').set(downLeftVertex3);
model.component('comp1').cpl('linext1').selection('srcvertex3').set(upperrightVertex4);
model.component('comp1').cpl('linext1').selection('dstvertex3').set(upperrightVertex3);


model.component('comp1').cpl('linext2').selection('srcvertex1').set(upperLeftVertex3);
model.component('comp1').cpl('linext2').selection('dstvertex1').set(upperLeftVertex4);
model.component('comp1').cpl('linext2').selection('srcvertex2').set(downLeftVertex3);
model.component('comp1').cpl('linext2').selection('dstvertex2').set(downLeftVertex4);
model.component('comp1').cpl('linext2').selection('srcvertex3').set(upperrightVertex3);
model.component('comp1').cpl('linext2').selection('dstvertex3').set(upperrightVertex4);


model.physics().create("ht", "HeatTransferInSolidsAndFluids", "geom1");
model.physics('ht').selection.set([domainIndices_liquid_tot(:)',domainIndices_fin_in_liquid(:)']);
model.physics('ht').feature('fluid1').selection.set(domainIndices_liquid_tot);
model.component("comp1").physics("ht").feature("solid1").set("k_mat", "userdef");
model.component("comp1").physics("ht").feature("solid1").set("k", ["k_s", 0, 0, 0, "k_s", 0, 0, 0, "k_s"]);
model.component("comp1").physics("ht").feature("solid1").set("rho_mat", "userdef");
model.component("comp1").physics("ht").feature("solid1").set("rho", "rho_s");
model.component("comp1").physics("ht").feature("solid1").set("Cp_mat", "userdef");
model.component("comp1").physics("ht").feature("solid1").set("Cp", "cp_s");
model.component("comp1").physics("ht").feature("fluid1").set("u_src", "userdef");
model.component("comp1").physics("ht").feature("fluid1").set("u", ["x0*u", "x0*v", "0"]);
model.component("comp1").physics("ht").feature("fluid1").set("fluidType", "gasLiquid");
model.component("comp1").physics("ht").feature("fluid1").set("rho_mat", "userdef");
model.component("comp1").physics("ht").feature("fluid1").set("rho", "rho_f");
model.component("comp1").physics("ht").feature("fluid1").set("Cp_mat", "userdef");
model.component("comp1").physics("ht").feature("fluid1").set("Cp", "cp_f");
model.component("comp1").physics("ht").feature("fluid1").set("k_mat", "userdef");
model.component("comp1").physics("ht").feature("fluid1").set("k", ["k_f", "0", "0", "0", "k_f", "0", "0", "0", "k_f"]);
model.component("comp1").physics("ht").create("hs1", "HeatSource", 2);
model.component("comp1").physics("ht").create("hs2", "HeatSource", 2);
model.physics('ht').feature('hs1').selection.set(domainIndices_liquid_solid);
model.physics('ht').feature('hs2').selection.set(domainIndices_fin_in_liquid);
model.component("comp1").physics("ht").feature("hs1").set("Q0", "x0*h_f*(comp1.linext1(T2)-T)/(dz_channel)");
model.component("comp1").physics("ht").feature("hs2").set("Q0", "x0*h_s*(comp1.linext1(T2)-T)/(dz_channel)");
model.component("comp1").physics("ht").create("ifl1", "Inflow", 1);
model.component("comp1").physics("ht").feature("ifl1").selection().set(1);
model.component("comp1").physics("ht").feature("ifl1").set("Tustr", "T_inlet");
model.component("comp1").physics("ht").create("ofl1", "ConvectiveOutflow", 1);
model.component("comp1").physics("ht").feature("ofl1").selection().set(14);
model.component("comp1").physics("ht").create("pc1", "PeriodicHeat", 1);
model.component("comp1").physics("ht").feature("pc1").selection().set([2, 3, 8, 9, 12 ,13]);


model.physics().create("spf", "LaminarFlow", "geom1");
model.physics("spf").selection.set(domainIndices_liquid_tot);
model.component("comp1").physics("spf").feature("fp1").set("rho_mat", "userdef");
model.component("comp1").physics("spf").feature("fp1").set("rho", "rho_f");
model.component("comp1").physics("spf").feature("fp1").set("mu_mat", "userdef");
model.component("comp1").physics("spf").feature("fp1").set("mu", "mu_f/x0");
model.component("comp1").physics("spf").feature("init1").set("u_init", ["u_inlet", "0", "0"]);
model.component("comp1").physics("spf").create("inl1", "InletBoundary", 1);
model.component("comp1").physics("spf").feature("inl1").selection().set(1);
model.component("comp1").physics("spf").feature("inl1").set("ComponentWise", "VelocityFieldComponentWise");
model.component("comp1").physics("spf").feature("inl1").set("u0", ["u_inlet" "0" "0"]);
model.component("comp1").physics("spf").create("out1", "OutletBoundary", 1);
model.component("comp1").physics("spf").feature("out1").selection().set(14);
model.component("comp1").physics("spf").create("sym1", "Symmetry", 1);
model.component("comp1").physics("spf").feature("sym1").selection().set([2, 3, 8, 9, 12, 13]);



model.physics().create("ht2", "HeatTransfer", "geom1");
model.physics('ht2').selection.set(domainIndices_down_rectangle);
model.component("comp1").physics("ht2").feature("solid1").set("Cp_mat", "userdef");
model.component("comp1").physics("ht2").feature("solid1").set("Cp", "cp_s");
model.component("comp1").physics("ht2").feature("solid1").set("rho_mat", "userdef");
model.component("comp1").physics("ht2").feature("solid1").set("rho", "rho_s");
model.component("comp1").physics("ht2").feature("solid1").set("k_mat", "userdef");
model.component("comp1").physics("ht2").feature("solid1").set("k", ["k_s", 0, 0, 0, "k_s", 0, 0, 0, "k_s"]);
model.component("comp1").physics("ht2").create("hs1", "HeatSource", 2);
model.component("comp1").physics("ht2").create("hs2", "HeatSource", 2);
model.component("comp1").physics("ht2").feature("hs1").selection().set(domainIndices_solid_liquid);
model.component("comp1").physics("ht2").feature("hs1").set("Q0", "x0*(Q0)/(dz_bp)-x0*h_f*(T2-comp1.linext2(T))/(dz_bp)");
model.component("comp1").physics("ht2").feature("hs2").selection().set(domainIndices_fin_in_solid);
model.component("comp1").physics("ht2").feature("hs2").set("Q0", "x0*(Q0)/(dz_bp)-x0*h_s*(T2-comp1.linext2(T))/(dz_bp)");


model.component("comp1").mesh().create("mesh1");
model.component("comp1").mesh("mesh1").feature("size").set("custom", false);
model.component("comp1").mesh("mesh1").feature("size").set("table", "cfd");
model.component("comp1").mesh("mesh1").feature("size").set("hauto", 1);
model.component("comp1").mesh("mesh1").create("ftri1", "FreeTri");
model.component("comp1").mesh("mesh1").feature("ftri1").selection().set([allDomains(:)',allDomains2(:)',domainIndices_right_rectangle(:)',...
    domainIndices_left_rectangle(:)',domainIndices_copy3(:)',domainIndices_copy4(:)']);
model.mesh('mesh1').run();






% Get the adjacency matrix (boundaries connected to domains)
boundaryIDs = mphgetadj(model, 'geom1', 'boundary', 'domain');

% Retrieve boundaries associated with the given domain ID
boundaries = boundaryIDs{domainIndices_left_rectangle}; % Replace 'domainID' with your specific domain
boundaries2 = boundaryIDs{domainIndices_right_rectangle};
% disp('Boundaries associated with the domain:');
% disp(boundaries2);
% disp(boundaries);

minX=Inf;
% Loop through the boundaries to find the one with the smallest x-coordinate
he_bound=zeros(2);j=1;
for i = 1:length(boundaries)
    % Get the coordinates of the current boundary using mphgetcoords
    coords = mphgetcoords(model, 'geom1', 'boundary', boundaries(i));
    if (coords(1,1)==coords(1,2))&&(coords(1,1)<minX)
        minX=coords(1,1);
        leftmostBoundary = boundaries(i);
    end
end

maxX=-Inf;
for i = 1:length(boundaries2)
    % Get the coordinates of the current boundary using mphgetcoords  
    coords2 = mphgetcoords(model, 'geom1', 'boundary', boundaries2(i));
    if (coords2(1,1)==coords2(1,2))&&(coords2(1,1)>maxX)
        maxX=coords2(1,1);
        rightmostBoundary = boundaries2(i);
    end
end


% Create an Integration coupling operator for the line
model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').selection.geom('geom1', 1); % 1 = boundary (line)
model.component('comp1').cpl('intop1').selection.set(leftmostBoundary); % Set the boundary ID

model.component('comp1').cpl.create('aveop1', 'Average');
model.component('comp1').cpl('aveop1').selection.geom('geom1', 1); % 1 = boundary (line)
model.component('comp1').cpl('aveop1').selection.set(leftmostBoundary); % Replace with actual boundary IDs


model.component('comp1').cpl.create('intop2', 'Integration');
model.component('comp1').cpl('intop2').selection.geom('geom1', 1); % 1 = boundary (line)
model.component('comp1').cpl('intop2').selection.set(rightmostBoundary); % Replace with actual boundary IDs

model.component('comp1').cpl.create('aveop2', 'Average');
model.component('comp1').cpl('aveop2').selection.geom('geom1', 2); 
model.component('comp1').cpl('aveop2').selection.set(domainIndices_copy3); % Replace with actual boundary IDs


model.component('comp1').cpl.create('aveop3', 'Average');
model.component('comp1').cpl('aveop3').selection.geom('geom1', 2); 
model.component('comp1').cpl('aveop3').selection.set(domainIndices_copy4); % Replace with actual boundary IDs


model.component('comp1').cpl.create('aveop4', 'Average');
model.component('comp1').cpl('aveop4').selection.geom('geom1', 2); 
model.component('comp1').cpl('aveop4').selection.set([domainIndices_fin_in_solid(:)',domainIndices_solid_liquid(:)']); % Replace with actual boundary IDs

model.component('comp1').cpl.create('aveop5', 'Average');
model.component('comp1').cpl('aveop5').selection.geom('geom1', 1); % 1 = boundary (line)
model.component('comp1').cpl('aveop5').selection.set(rightmostBoundary); % Replace with actual boundary IDs

model.component('comp1').cpl.create('maxop1', 'Maximum');
model.component('comp1').cpl('maxop1').selection.geom('geom1', 2); % 2 = area (2-dimensional)
model.component('comp1').cpl('maxop1').selection.set([domainIndices_fin_in_solid(:)',domainIndices_solid_liquid(:)']); 


model.component('comp1').cpl.create('maxop4', 'Maximum');
model.component('comp1').cpl('maxop4').selection.geom('geom1', 2); % 2 = area (2-dimensional)
model.component('comp1').cpl('maxop4').selection.set(domainIndices_copy3);



model.component('comp1').cpl.create('maxop5', 'Maximum');
model.component('comp1').cpl('maxop5').selection.geom('geom1', 2); % 2 = area (2-dimensional)
model.component('comp1').cpl('maxop5').selection.set(domainIndices_copy4);



model.variable.create('varGlobal');
model.variable('varGlobal').label('GlobalVariables');
model.variable('varGlobal').set('mass_in', 'rho_f*H*dz_channel*u_inlet');
model.variable('varGlobal').set('T_out', 'Q0*H^2/(mass_in*cp_f)+T_inlet');
model.variable('varGlobal').set('T_s_ave','comp1.aveop4(T2)');
model.variable('varGlobal').set('T_f_ave', '(T_out+T_inlet)/2');
model.variable('varGlobal').set('x_center1', 'comp1.aveop2(x)');
model.variable('varGlobal').set('y_center1', 'comp1.aveop2(y)');
model.variable('varGlobal').set('x_center2', 'comp1.aveop3(x)');
model.variable('varGlobal').set('y_center2', 'comp1.aveop3(y)');
model.variable('varGlobal').set('T_out_2', 'comp1.intop2(u*T)/comp1.intop2(u)');
model.variable('varGlobal').set('T_max', 'comp1.maxop1(T2)');
model.variable('varGlobal').set('pressure_drop', 'comp1.aveop5(p)-comp1.aveop1(p)');

model.variable('varGlobal').set('max_x_fin1','comp1.maxop4(abs(x - x_center1))');
model.variable('varGlobal').set('max_y_fin1','comp1.maxop4(abs(y - y_center1))');
model.variable('varGlobal').set('max_x_fin2','comp1.maxop5(abs(x - x_center2))');
model.variable('varGlobal').set('max_y_fin2','comp1.maxop5(abs(y - y_center2))');



geom.run();

model.study().create("std1");
model.study("std1").create("stat", "Stationary");
model.study('std1').run();


T_max=mphglobal(model, 'T_max');
%disp('T_max is');disp(T_max);
Ts_ave=mphglobal(model, 'T_s_ave');
%disp('Ts_ave is');disp(Ts_ave);
pressure_drop=mphglobal(model, 'pressure_drop');
%disp('pressure_drop is');disp(pressure_drop);
T_out=mphglobal(model, 'T_out');
%disp('T_out is');disp(T_out);


T_out_2=mphglobal(model, 'T_out_2');
%disp('T_out_2 is');disp(T_out_2);

% Initialize arrays for x_cen_fin and y_cen_fin
x_cen_fin = zeros(num_fins, 1);
y_cen_fin = zeros(num_fins, 1);

% Evaluate center coordinates
x_cen_fin(1) = mphglobal(model, 'x_center1');
x_cen_fin(2) = mphglobal(model, 'x_center2');

y_cen_fin(1) = mphglobal(model, 'y_center1');
y_cen_fin(2) = mphglobal(model, 'y_center2');


max_x_fin1=mphglobal(model,'max_x_fin1');
max_x_fin2=mphglobal(model,'max_x_fin2');
max_y_fin1=mphglobal(model,'max_y_fin1');
max_y_fin2=mphglobal(model,'max_y_fin2');
max_deformation = zeros(num_fins, 1);

% Compute max_deformation
max_deformation(1) = min(max_x_fin1, max_y_fin1);
max_deformation(2) = min(max_x_fin2, max_y_fin2);
% max_deformation(3) = min(max_x_fin3, max_y_fin3); % Uncomment if using max_x_fin3 and max_y_fin3

% Initialize arrays for x_range and y_range
x_range = zeros(num_fins, 1);
y_range = zeros(num_fins, 1);

% Compute x_range and y_range
%x_range(1) = (fluid_solid_area_width / 2) - 1.05 * max_x_fin1;
%x_range(2) = (fluid_solid_area_width / 2) - 1.05 * max_x_fin2;
max_deform=0.5;
x_range(1) = (fluid_solid_area_width / 2) - 1.05 * max_deform;
x_range(2) = (fluid_solid_area_width / 2) - 1.05 * max_deform;
% x_range(3) = (fluid_solid_area_width / 2) - 1.05 * max_x_fin3; % Uncomment if using max_x_fin3

%y_range(1) = (fluid_solid_area_height / 2) - 1.05 * max_y_fin1;
%y_range(2) = (fluid_solid_area_height / 2) - 1.05 * max_y_fin2;
y_range(1) = (fluid_solid_area_height / 2) - 1.05 * max_deform;
y_range(2) = (fluid_solid_area_height / 2) - 1.05 * max_deform;
% y_range(3) = (fluid_solid_area_height / 2) - 1.05 * max_y_fin3; % Uncomment if using max_y_fin3


objective_both = ((Ts_ave - T_inlet) + (T_max - T_inlet));
objective_ave = (Ts_ave - T_inlet);
objective_max = (T_max- T_inlet);
%disp('objective_ave is');disp(objective_ave);


model.result().create("pg1", "PlotGroup2D");
model.result('pg1').set('data', 'dset1');
model.result("pg1").create("surf1", "Surface");
model.result("pg1").feature("surf1").set("expr", "T");
model.result("pg1").create("surf2", "Surface");
model.result("pg1").feature("surf2").set("expr", "T2");





model.result().export().create("plot1", "Plot");
model.result().export("plot1").set("plot", "surf1");
model.result().export("plot1").set("exporttype", "vtu");
filename = sprintf('/home/morteza/Desktop/New Folder 6/pmcopy3/T_distribution_%d.vtu', iteration);
%filename = sprintf('/home/morteza/Desktop/T_distribution/T_distribution_1.vtu');
model.result().export("plot1").set("filename", filename);
model.result().export("plot1").run();

model.result().export().create("plot2", "Plot");
model.result().export("plot2").set("plot", "surf2");
model.result().export("plot2").set("exporttype", "vtu");
filename = sprintf('/home/morteza/Desktop/New Folder 6/pmcopy3/T2_distribution_%d.vtu', iteration);

%filename = sprintf('/home/morteza/Desktop/T_distribution/T2_distribution_1.vtu');
model.result().export("plot2").set("filename", filename);
model.result().export("plot2").run();





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model.result().export().create("data1", "Data");
% model.result().export("data1").remove("unit", 0);
% model.result().export("data1").remove("descr", 0);
% model.result().export("data1").remove("expr", new int[]{0});
% model.result().export("data1").remove("unit", 0);
% model.result().export("data1").remove("descr", 0);
% model.result().export("data1").remove("expr", new int[]{0});
% model.result().export("data1").setIndex("expr", 1, 0);
% model.result().export("data1").set("data", "edg1");
% model.result().export("data1").set("filename", "/home/morteza/Desktop/Untitled.txt");
% 
% model.result().dataset("edg1").selection().set(17, 18, 22, 25);
% model.result().dataset("edg1").selection().named("geom1_cc1_bnd");
% model.result().dataset("edg1").selection().named("geom1_cc2_bnd");
% model.result().dataset("edg1").selection().set();
% model.result().dataset("edg1").selection().named("geom1_cc2_bnd");
% 
% model.result().dataset().create("edg1", "Edge2D");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model.result().dataset().create("edg1", "Edge2D");
model.result().dataset("edg1").set("data", "dset1");
model.result().dataset("edg1").selection().named("geom1_cc1_bnd");

model.result().dataset().create("edg2", "Edge2D");
model.result().dataset("edg2").set("data", "dset1");
model.result().dataset("edg2").selection().named("geom1_cc2_bnd");

model.result().dataset().create("edg3", "Edge2D");
model.result().dataset("edg3").set("data", "dset1");
model.result().dataset("edg3").selection().named("geom1_fin_in_liquid_bnd");


% model.result().export("plot2").set("plot", "surf2");
% model.result().export("plot2").set("exporttype", "vtu");
% model.result().export("plot2").set("filename", "/home/morteza/Desktop/T_distribution/T2_distribution_3.vtu");



% model.result().export("data1").remove("unit", 0);
% model.result().export("data1").remove("descr", 0);
% 
% model.result().export("data1").remove("unit", 0);
% model.result().export("data1").remove("descr", 0);
% 
% model.result().export("data1").setIndex("expr", 1, 0);
% model.result().export("data1").set("data", "edg1");
% model.result().export("data1").set("filename", "/home/morteza/Desktop/Untitled.txt");

model.result().export().create("Union", "Data");
model.result().export("Union").setIndex("expr", 1, 0);
model.result().export("Union").setIndex("expr", "", 1);
model.result().export("Union").set("data", "edg3");
filename = sprintf('/home/morteza/Desktop/New Folder 6/pmcopy3/Union_%d.plt', iteration);
model.result().export("Union").set("filename", filename);
model.result().export("Union").run();


model.result().export().create("bez_curve_1", "Data");
model.result().export("bez_curve_1").setIndex("expr", 1, 0);
model.result().export("bez_curve_1").setIndex("expr", "", 1);
model.result().export("bez_curve_1").set("data", "edg1");
filename = sprintf('/home/morteza/Desktop/New Folder 6/pmcopy3/bez1_%d.plt', iteration);
model.result().export("bez_curve_1").set("filename", filename);
model.result().export("bez_curve_1").run();


model.result().export().create("bez_curve_2", "Data");
model.result().export("bez_curve_2").setIndex("expr", 1, 0);
model.result().export("bez_curve_2").setIndex("expr", "", 1);
model.result().export("bez_curve_2").set("data", "edg2");
filename = sprintf('/home/morteza/Desktop/New Folder 6/pmcopy3/bez2_%d.plt', iteration);
model.result().export("bez_curve_2").set("filename", filename);
model.result().export("bez_curve_2").run();


%savePath = 'pmcopy3_two_shape_uniform_Heat_source_2.mph';
%model.save(savePath);

% --- after you have pulled Ts_ave, T_max, etc. ---

% 1) throw away *all* solution vectors --------------------- 👉
model.sol('sol1').clearSolutionData;          % 5.4+ API
% or, in older versions: model.sol.remove('sol1');

% 2) wipe cached Results data ------------------------------ 👉
model.result().clearStoredPlotData();         % 5.5+  (docs ref)




for t = ModelUtil.tags
    ModelUtil.remove(t);           % deregister
end
clearvars -regexp ^model           % clear model1, model2, …
java.lang.System.gc();

ModelUtil.remove('Model1');
clear model  
java.lang.System.gc();



