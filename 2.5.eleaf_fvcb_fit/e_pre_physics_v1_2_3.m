% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.5

load save_e_geom.mat %from e_geo_main_v1_2.m
copyfile('../2.e_raytracing/results_cell*','./');

%% readin results from ray tracing
%tmp_filename = input('Ray tracing results (e.g. results_cell):','s');
%if isempty(tmp_filename)
    tmp_filename = 'results_cell_470';
%end
file_ID=fopen(tmp_filename,'r');
load save_raytracing.mat count_mschl%from geo_export_e_geo_main_v1_2.m
formatSpec = '%d %d %d %d %d %e';
cutnum_x=500;cutnum_y=25;tmp_ab_profile=zeros(2,size(count_mschl,2));
for loop_i=1:size(count_mschl,2)
    [tmp_read,count]=fscanf(file_ID,formatSpec,6);%cell wall
    [tmp_read,count]=fscanf(file_ID,formatSpec,6*count_mschl(loop_i));%chl
    tmp_read=reshape(tmp_read,6,[]);
    tmp_ab_profile(1,loop_i)=sum(tmp_read(6,:))/(cutnum_x*cutnum_y);
    [tmp_read,count]=fscanf(file_ID,formatSpec,6);%vac
end
fclose(file_ID);
%ab_all=sum(tmp_ab_profile(1,:));


%tmp_filename = input('Ray tracing results (e.g. results_cell):','s');
%if isempty(tmp_filename)
    tmp_filename = 'results_cell_665';
%end
file_ID=fopen(tmp_filename,'r');
load save_raytracing.mat count_mschl%from geo_export_e_geo_main_v1_2.m
formatSpec = '%d %d %d %d %d %e';
cutnum_x=500;cutnum_y=25;%tmp_ab_profile=zeros(2,size(count_mschl));
for loop_i=1:size(count_mschl,2)
    [tmp_read,count]=fscanf(file_ID,formatSpec,6);%cell wall
    [tmp_read,count]=fscanf(file_ID,formatSpec,6*count_mschl(loop_i));%chl
    tmp_read=reshape(tmp_read,6,[]);
    tmp_ab_profile(2,loop_i)=sum(tmp_read(6,:))/(cutnum_x*cutnum_y);
    [tmp_read,count]=fscanf(file_ID,formatSpec,6);%vac
end
fclose(file_ID);
ab_all=sum(tmp_ab_profile(2,:),2);

ab_profile=[0.1,0.9]*tmp_ab_profile;
bar(tmp_ab_profile')
save('ab_profile.mat','ab_profile')
