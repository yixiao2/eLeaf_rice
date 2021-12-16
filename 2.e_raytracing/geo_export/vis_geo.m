% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.4

clear;clc
colormap([0.3,0.3,0.3;0.2,0.8,0.2;0.2,0.2,0.8]);
alpha=0.2;
linecolor=[0.8,0.8,0.8];
wall=0;
chl=1;
vac=2;

load save_raytracing count_nscell
for loop_i=1:4
    tmp_name=['ns',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read_xy(tmp_name);
    trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
end

load save_raytracing count_mschl count_MSC
for loop_i=1:count_MSC
    %if loop_i==1
    tmp_name=['ms',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read_xy(tmp_name);
    trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    for loop_j=1:count_mschl(loop_i)
        tmp_name=['ms',num2str(loop_i),'c',num2str(loop_j),'.ply'];
        [tri,pts]=ply_read_xy(tmp_name);
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*chl,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    end
    tmp_name=['ms',num2str(loop_i),'v.ply'];
    [tri,pts]=ply_read_xy(tmp_name);
    trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*vac,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    %end
end
axis equal;

%% plot RT_file4plto
tmp_paths_all=load('../test_plot4RT');%%need to change Defs.h to get this file
tmp_idx=find(tmp_paths_all(:,1)==80&tmp_paths_all(:,2)==10);
plot3(tmp_paths_all(tmp_idx,3),tmp_paths_all(tmp_idx,4),tmp_paths_all(tmp_idx,5),'r-o','linewidth',1);hold on;
axis equal
view(90,0)

% P=[3.53688e-005,1.62161e-005,0];
% D=[0.67129415129712866, -0.52763743027764065, -0.52054097303120694];
% dis=1e-5;
% plot3([P(1);P(1)+D(1)*dis],[P(2);P(2)+D(2)*dis],[P(3);P(3)+D(3)*dis],'r-o','linewidth',2);hold on;

% P=[4.63498e-005,1.92781e-005,0];
% D=[4.62101e-005,1.93588e-005,0];
% dis=1e-15;
% plot3([P(1);D(1)],[P(2);D(2)],[P(3);D(3)],'r-o','linewidth',2);hold on;

