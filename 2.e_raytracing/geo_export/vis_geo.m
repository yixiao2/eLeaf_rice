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

% %% read debug_RT
% ray_i_start=0;
% ray_i_end=500;
% ray_j_start=0;
% ray_j_end=25;
% [fid,Msg]=fopen('debug_RT.txt','rt');
% for i=ray_i_start:ray_i_end-1
%     for j=ray_j_start:ray_j_end-1
%         tmp_num_pts=fscanf(fid,'%d\n',1);
%         clear tmp_tri
%         for k=1:tmp_num_pts
%             tmp_tri(k,:)=fscanf(fid,'%e %e %e\n',3);
%         end
%         fscanf(fid,'\n');
%         %if(i>90&&i<92)
%         plot3(tmp_tri(:,1),tmp_tri(:,2),tmp_tri(:,3),'r-o','linewidth',2);hold on;
%         %end
%     end
% end
% fclose(fid);



% P=[3.53688e-005,1.62161e-005,0];
% D=[0.67129415129712866, -0.52763743027764065, -0.52054097303120694];
% dis=1e-5;
% plot3([P(1);P(1)+D(1)*dis],[P(2);P(2)+D(2)*dis],[P(3);P(3)+D(3)*dis],'r-o','linewidth',2);hold on;

P=[4.63498e-005,1.92781e-005,0];
D=[4.62101e-005,1.93588e-005,0];
dis=1e-15;
plot3([P(1);D(1)],[P(2);D(2)],[P(3);D(3)],'r-o','linewidth',2);hold on;

