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


model =mphload('pmcopy3_two_shape_uniform_Heat_source.mph');


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



model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "x11", 0, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "y11", 1, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "x12", 0, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "y12", 1, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "x13", 0, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "y13", 1, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "x14", 0, 3);
model.component("comp1").geom("geom1").feature("cc1").feature("cb1").setIndex("p", "y14", 1, 3);



model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "x14", 0, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "y14", 1, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "x15", 0, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "y15", 1, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "x16", 0, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "y16", 1, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "x17", 0, 3);
model.component("comp1").geom("geom1").feature("cc1").feature("cb2").setIndex("p", "y17", 1, 3);


model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "x17", 0, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "y17", 1, 0);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "y18", 1, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "x18", 0, 1);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "x19", 0, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "y19", 1, 2);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "x110", 0, 3);
model.component("comp1").geom("geom1").feature("cc1").feature("cb3").setIndex("p", "y110", 1, 3);


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



model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "x21", 0, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "y21", 1, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "x22", 0, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "y22", 1, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "x23", 0, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "y23", 1, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "x24", 0, 3);
model.component("comp1").geom("geom1").feature("cc2").feature("cb1").setIndex("p", "y24", 1, 3);


model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "x24", 0, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "y24", 1, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "x25", 0, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "y25", 1, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "x26", 0, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "y26", 1, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "x27", 0, 3);
model.component("comp1").geom("geom1").feature("cc2").feature("cb2").setIndex("p", "y27", 1, 3);


model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "x27", 0, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "y27", 1, 0);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "y28", 1, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "x28", 0, 1);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "x29", 0, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "y29", 1, 2);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "x210", 0, 3);
model.component("comp1").geom("geom1").feature("cc2").feature("cb3").setIndex("p", "y210", 1, 3);


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


%model.result('pg1').set('data', 'dset1');
%model.result("pg1").feature("surf1").set("expr", "T");
%model.result("pg1").feature("surf2").set("expr", "T2");










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
