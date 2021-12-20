% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <xiaoyi@sippe.ac.cn>
% @version: 1.2.6

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

%% run geo_export to generate triangle meshed surface
cd ../2.e_raytracing/geo_export
geo_export_e_geo_main_v1_2
%generate Defs.h in geo_export
copyfile('Defs_template.h','Defs.h');
fileID=fopen('Defs.h','a');
load save_e_geom.mat xmax xmin ymax ymin zmax zmin
fprintf(fileID,'#define xmax %e\n#define xmin %e\n#define zmax %e\n#define zmin %e\n#define ymax %e\n#define ymin %e\n\n',xmax,xmin,zmax,zmin,ymax,ymin);
fprintf(fileID,'#define bsc_num 1\n#define msc_num %d\n#define msall_num (bsc_num+msc_num)\n',count_MSC);
fprintf(fileID,'int count_ms;\nint count_nonms;\n#define ms_max_chl_num %d\nint ms_chl_num[msall_num];\ndouble ms_chl_con[msall_num];\n#define nonms_num %d\n\n',max(count_mschl),count_nscell);
fprintf(fileID,'Object *p_cell_ms[msall_num];\nObject *p_chl_ms[msall_num][ms_max_chl_num];\nObject *p_vac_ms[msall_num];\nObject *p_leaf;\nObject *p_cell_ns[nonms_num];\n\n');
fprintf(fileID,'#endif\n');
fclose(fileID);
copyfile('Defs.h','../');
copyfile('count_chl4RT','../');%%file generated from geo_export.m

%copy ms*.ply to ../MS/
copyfile('ms*.ply','../MS/');
%copy ns*.ply to ../nonMS/
copyfile('ns*.ply','../nonMS/');

copyfile('save_raytracing.mat','../../2.5.eleaf_fvcb_fit');

%% ray tracing
cd ..
system('make');
%system('sh run_raytracing_475.sh');
%system('sh run_raytracing_625.sh');
SAC_licor_blue475nm=[1.14e-4*100, 4.26e4*1e-4];
SAC_licor_red625nm=[2.834e-3*100, 2.34e4*1e-4];
tmp_SAC_water=min(SAC_licor_blue475nm(1),SAC_licor_red625nm(1));
tmp_SAC_chl=min(SAC_licor_blue475nm(2),SAC_licor_red625nm(2));
%%%%%%%%% CASE5: cut x&y every 5 rays; batch number 100*100
rep_num=1;
CASE_NUM=5;
%delete(gcp('nocreate'));
%parpool(1);% num of cores for ray tracing
RT_x_range=500;
RT_y_range=25;
RT_x_perthread=5;
RT_y_perthread=1;
num_loop_x=RT_x_range/RT_x_perthread;
num_loop_y=RT_y_range/RT_y_perthread;
count_batch=0;
batch_cmd_all={};
for loop_x=1:num_loop_x
    for loop_y=1:num_loop_y
        count_batch=count_batch+1;
        ray_x_start=RT_x_perthread*(loop_x-1);
        ray_x_end=RT_x_perthread*(loop_x);
        ray_y_start=RT_y_perthread*(loop_y-1);
        ray_y_end=RT_y_perthread*(loop_y);
        %%%% For fitting, run RT using 1/20*chl_SAC, assuming 1/10 is the lower limit of
        %%%% [chl] and 1/2 is the lower limit of chl_SAC itself.
        %%%% 1.22e4*1e-4/20=0.0610
        tmp_bash_cmd=['./test ',num2str(tmp_SAC_water),' ',num2str(tmp_SAC_chl),' ',...
            'results_abevents_tmpnm_500x_rep',num2str(rep_num),'_',num2str(count_batch),' ',...
            'count_chl4RT ',...
            num2str(RT_x_range),' ',num2str(RT_y_range),' ',...
            num2str(ray_x_start),' ',num2str(ray_x_end),' ',...
            num2str(ray_y_start),' ',num2str(ray_y_end),' ',...
            'results_srf_tmpnm_500x_rep',num2str(rep_num),'_',num2str(count_batch),' ',...
            'results_sum_tmpnm_500x_rep',num2str(rep_num),'_',num2str(count_batch),...
            ' >>results_RTlog_tmpnm_500x_',num2str(rep_num),'_',num2str(count_batch)];
        %system(['start ',tmp_batch_cmd])
        batch_cmd_all={batch_cmd_all{:},tmp_bash_cmd};
    end
end
tic
%parfor loop_bash=1:size(batch_cmd_all,2)
for loop_bash=1:size(batch_cmd_all,2)
    system(batch_cmd_all{loop_bash});
end
time_case(rep_num)=toc;
%%%% check results_sum files
count_batch=0;
record_success=zeros(1,num_loop_x*num_loop_y);
for loop_x=1:num_loop_x
    for loop_y=1:num_loop_y
        count_batch=count_batch+1;
        tmp_file_name=['results_sum_tmpnm_500x_rep',num2str(rep_num),'_',num2str(count_batch)];
        if(exist(tmp_file_name,'file')==2)
            record_success(count_batch)=1;
        end
    end
end
failed_threads{rep_num}=find(record_success==0);

if numel(failed_threads{rep_num})/(num_loop_x*num_loop_y)>0.05
    %%error('High fail rate during ray tracing.')
    disp('[Warning]: High fail rate during ray tracing.')
end

%% trace_recal --> ab profiles under 475nm and 625nm
num_layer4lightprofile=10;
%% SAC_licor_blue475nm=[1.14e-4*100, 4.26e4*1e-4];
RT_x_range=500;
RT_y_range=25;
RT_x_perthread=5;
RT_y_perthread=1;
num_loop_x=RT_x_range/RT_x_perthread;
num_loop_y=RT_y_range/RT_y_perthread;
tmp_bash_cmd=['./trace_recal ',num2str(SAC_licor_blue475nm(1)),' ',num2str(SAC_licor_blue475nm(2)),' count_chl4RT ',...
    num2str(RT_x_range),' ',num2str(RT_y_range),' ',...
    num2str(RT_x_perthread),' ',num2str(RT_y_perthread),' ',...
    'results_abevents_tmpnm_500x_rep',num2str(rep_num),'_ ',...
    'results_sum_tmpnm_500x_rep',num2str(rep_num),'_ ',...
    num2str(num_layer4lightprofile),' ',...
    'results_merged_abtri_475nm_500x_rep',num2str(rep_num),' ',...
    'results_merged_absrf_475nm_500x_rep',num2str(rep_num),' ',...
    'results_merged_abprofile_475nm_500x_rep',num2str(rep_num),'_layerN',num2str(num_layer4lightprofile),' ',...
    'results_merged_rtsum_475nm_500x_rep',num2str(rep_num)];
system(tmp_bash_cmd);
%% SAC_licor_red625nm=[2.834e-3*100, 2.34e4*1e-4];
RT_x_range=500;
RT_y_range=25;
RT_x_perthread=5;
RT_y_perthread=1;
num_loop_x=RT_x_range/RT_x_perthread;
num_loop_y=RT_y_range/RT_y_perthread;
tmp_bash_cmd=['./trace_recal ',num2str(SAC_licor_red625nm(1)),' ',num2str(SAC_licor_red625nm(2)),' count_chl4RT ',...
    num2str(RT_x_range),' ',num2str(RT_y_range),' ',...
    num2str(RT_x_perthread),' ',num2str(RT_y_perthread),' ',...
    'results_abevents_tmpnm_500x_rep',num2str(rep_num),'_ ',...
    'results_sum_tmpnm_500x_rep',num2str(rep_num),'_ ',...
    num2str(num_layer4lightprofile),' ',...
    'results_merged_abtri_625nm_500x_rep',num2str(rep_num),' ',...
    'results_merged_absrf_625nm_500x_rep',num2str(rep_num),' ',...
    'results_merged_abprofile_625nm_500x_rep',num2str(rep_num),'_layerN',num2str(num_layer4lightprofile),' ',...
    'results_merged_rtsum_625nm_500x_rep',num2str(rep_num)];
system(tmp_bash_cmd);

%% tar intermediate files from ray tracing
system('tar -zcf results_files_abevents.tar.gz results_abevents_* --remove-files');% tar -zxvf ***.tar.gz; -v = output log
system('tar -zcf results_files_srf.tar.gz results_srf_* --remove-files');
system('tar -zcf results_files_sum.tar.gz results_sum_* --remove-files');
system('tar -zcf results_files_RTlog.tar.gz results_RTlog_* --remove-files');

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
