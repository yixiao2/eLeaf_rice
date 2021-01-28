% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <xiaoyi@sippe.ac.cn>
% @version: 1.2.5

function run_e_leaf_v1_2_5(CFG_MODE,CFG_PARA_COM,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - CFG_MODE=1: before fit
% - CFG_MODE=2: dissection
% - CFG_MODE=3: sensitivity analysis
% - CFG_MODE=4: after fit; for std of fit
% - CFG_PARA_COM= e.g. [1 1 0 0 0 0 0 0 0]
% optional CFG_SA_SET,CFG_SA_FC
% - CFG_SA_SET= variable names
% - CFG_SA_FC= scale factor for CFG_SA_SET
% - run_e_leaf_v1_2_5 is a function with parameters above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

%%%%%%%%%%%%open comsol server port%%%%%%%%%%
%N_pool=12;
%tmp_port=2035+[1:13,15:25];
%for tmp_loop=1:N_pool
%    comsolPort=tmp_port(tmp_loop);
%    system( ['gnome-terminal -x comsol -np 2 server -silent -port ' num2str(comsolPort)] );
%end

%% parse parameters
p = inputParser;
p.StructExpand=true; % if we allow parameter  structure expanding
p.addRequired('CFG_MODE', @isnumeric);
p.addRequired('CFG_PARA_COM', @isnumeric);
p.addOptional('CFG_SA_SET', @iscellstr);
p.addOptional('CFG_SA_FC', @isnumeric);
p.parse(CFG_MODE,CFG_PARA_COM,varargin{:});
r=p.Results;
CFG_SA_SET=r.CFG_SA_SET;
CFG_SA_FC=r.CFG_SA_FC;

save CFG.mat CFG*

if CFG_MODE==1
    e_geo_parainput_v1_2_5_b4fit(CFG_PARA_COM);
else
    e_geo_parainput_v1_2_5_a4tfit(CFG_PARA_COM);
    if CFG_MODE==3
        load parainput.mat
        for tmp_loop=1:length(CFG_SA_SET)
            tmp_chr=CFG_SA_SET{tmp_loop};
            if exist(tmp_chr,'var')==0
                error('SA error:no this variable - %s',tmp_chr)
            end
            eval([tmp_chr,'=',tmp_chr,'*',num2str(CFG_SA_FC(tmp_loop)),';'])
        end
        recal_parainput;
        save parainput.mat
    end
end

e_geo_MS3D_v1_2_1 % for 3D mesophyll cell
e_geo_IAS2D_v1_2_1 % for 2D distribution
e_geo_main_v1_2_1 % for 3D leaf geometry
e_geo_test_v1_2_3 % test geometry compatibility
e_physics_addcresel_v1_1 % identify domain and boundaries
AVR_cal % calculate volume of chloroplast

%% module of ray tracing
if(exist('../2.e_raytracing','dir')~=7)
    display('no dir 2.e_raytracing');
else
    copyfile('parainput.mat','../2.e_raytracing/geo_export');
    copyfile('tmp_IAS2D.mat','../2.e_raytracing/geo_export');
    copyfile('tmp_MS3D.mat','../2.e_raytracing/geo_export');
    copyfile('tmpCK_geomIP_cad_mesh_cresel.mph','../2.e_raytracing/geo_export');
    copyfile('tmpCK_geomIP_cad_mesh_cresel.mph','../2.5.eleaf_fvcb_fit');
    copyfile('save_e_geom.mat','../2.e_raytracing/geo_export');
    copyfile('save_e_geom.mat','../2.5.eleaf_fvcb_fit');
end

%generate leaf box (boundary of ray tracing)
fileID = fopen('leaf','w');
fprintf(fileID,'8\n12\n');
fprintf(fileID,'%e %e %e\n',xmin,ymax,zmin);
fprintf(fileID,'%e %e %e\n',xmax,ymax,zmin);
fprintf(fileID,'%e %e %e\n',xmax,ymin,zmin);
fprintf(fileID,'%e %e %e\n',xmin,ymin,zmin);
fprintf(fileID,'%e %e %e\n',xmin,ymax,zmax);
fprintf(fileID,'%e %e %e\n',xmax,ymax,zmax);
fprintf(fileID,'%e %e %e\n',xmax,ymin,zmax);
fprintf(fileID,'%e %e %e\n',xmin,ymin,zmax);
fprintf(fileID,'3 0 1 3\n3 1 2 3\n3 0 3 4\n3 3 4 7\n3 0 1 4\n3 1 4 5\n3 1 2 5\n3 2 5 6\n3 2 3 7\n3 2 6 7\n3 4 5 7\n3 5 6 7\n');
fclose(fileID);
copyfile('leaf','../2.e_raytracing');

%generate bash command; 470nm
fileID=fopen('run_raytracing_470.sh','w');
fprintf(fileID,'for i in {1..5}\ndo\n  if [ ! -f results_cell_470 ];\n  then\n');
fprintf(fileID,'./test 0 500 0 25 1.060000e-04 1.113930e+05 %e results_tri_470 results_cell_470 results_sum_470\n',chl_con4raytracing);
fprintf(fileID,'  else\n    break\n  fi\ndone\n\necho $i\nif [ ! -f results_cell_470 ];\n  then\n  echo failed\nfi\n');
fclose(fileID);
copyfile('run_raytracing_470.sh','../2.e_raytracing');

%generate bash command; 665nm
fileID=fopen('run_raytracing_665.sh','w');
fprintf(fileID,'for i in {1..5}\ndo\n  if [ ! -f results_cell_665 ];\n  then\n');
fprintf(fileID,'./test 0 500 0 25 4.290000e-03 5.277384e+04 %e results_tri_675 results_cell_665 results_sum_665\n',chl_con4raytracing);
fprintf(fileID,'  else\n    break\n  fi\ndone\n\necho $i\nif [ ! -f results_cell_665 ];\n  then\n  echo failed\nfi\n');
fclose(fileID);
copyfile('run_raytracing_665.sh','../2.e_raytracing');

%% run geo_export to generate triangle meshed surface
cd ../2.e_raytracing/geo_export
geo_export_e_geo_main_v1_2
%generate Defs.h in geo_export
copyfile('Defs_template.h','Defs.h');
fileID=fopen('Defs.h','a');
load save_e_geom.mat xmax xmin ymax ymin zmax zmin
fprintf(fileID,'#define xmax %e\n#define xmin %e\n#define zmax %e\n#define zmin %e\n#define ymax %e\n#define ymin %e\n\n',xmax,xmin,zmax,zmin,ymax,ymin);
fprintf(fileID,'#define ms_num %d\nint count_ms;\n#define ms_max_chl_num %d\nint ms_chl_num[ms_num];\n#define num_nonMS %d\n\n',count_MSC,max(count_mschl),count_nscell);
fprintf(fileID,'Object *p_cell_ms[ms_num];\nObject *p_chl_ms[ms_num][ms_max_chl_num];\nObject *p_vac_ms[ms_num];\nObject *p_leaf;\nObject *p_cell_ns[num_nonMS];\n\n');
fprintf(fileID,'#endif\n');
fclose(fileID);
copyfile('Defs.h','../');
copyfile('count_chl','../');

%copy ms*.ply to ../MS/
copyfile('ms*.ply','../MS/');
%copy ns*.ply to ../nonMS/
copyfile('ns*.ply','../nonMS/');

copyfile('save_raytracing.mat','../../2.5.eleaf_fvcb_fit');

%% ray tracing
cd ..
system('make');
system('sh run_raytracing_470.sh');
system('sh run_raytracing_665.sh');

%% run e-physics
load ../1.e_geom/CFG.mat
cd ../2.5.eleaf_fvcb_fit/
if CFG_MODE==1||CFG_MODE==4
    if all(CFG_PARA_COM==0)% CK
        ARG_TYPE=0;
        run_e_leaf_fvcb_fit;
    elseif all(CFG_PARA_COM==1)%HC
        ARG_TYPE=1;
        run_e_leaf_fvcb_fit;
    else
        error('Error in run_eleaf_fvcb_fit.m');
    end
elseif CFG_MODE==2
    run_e_leaf_fvcb_dissection
elseif CFG_MODE==3
    run_e_leaf_fvcb_SA
end

toc
