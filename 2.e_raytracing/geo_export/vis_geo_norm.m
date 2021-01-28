% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.4

clear;clc
colormap([0.3,0.3,0.3;0.2,0.8,0.2;0.2,0.2,0.8]);
alpha=1;
linecolor=[0.5,0.5,0.5];
wall=0;
chl=1;
vac=2;

% [tri,pts]=ply_read_xy('ep1.ply');
% trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
% 
% [tri,pts]=ply_read_xy('ms1.ply');
% trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
% [tri,pts]=ply_read_xy('ms1c1.ply');
% trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*chl,'facealpha',alpha,'EdgeColor',linecolor);hold on;
% [tri,pts]=ply_read_xy('ms1v.ply');
% trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*vac,'facealpha',alpha,'EdgeColor',linecolor);hold on;
% 
% [tri,pts]=ply_read_xy('ms2.ply');
% trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
% [tri,pts]=ply_read_xy('ms2c1.ply');
% trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*chl,'facealpha',alpha,'EdgeColor',linecolor);hold on;
% [tri,pts]=ply_read_xy('ms2v.ply');
% trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*vac,'facealpha',alpha,'EdgeColor',linecolor);hold on;
[tri,pts]=ply_read_xy('ms2c1.ply');
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),ones(size(pts,1),1)*wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
axis equal;

%% read debug_norvec
[fid,Msg]=fopen('debug_norvec.txt','rt');
vec_size=0.5;
for i=1:size(tri,1)
    tmp_o(i,:)=mean(pts(tri(i,:),:));
    tmp_vec(i,:)=fscanf(fid,'%e %e %e\n',3);
    %tmp_e(i,:)=tmp_o+tmp_vec'*vec_size;
    %plot3([tmp_o(1),tmp_e(1)],[tmp_o(2),tmp_e(2)],[tmp_o(3),tmp_e(3)],'r-');hold on;
    %quiver3(tmp_o(1),tmp_o(2),tmp_o(3),tmp_vec(1),tmp_vec(2),tmp_vec(3))
end
fclose(fid);
quiver3(tmp_o(:,1),tmp_o(:,2),tmp_o(:,3),tmp_vec(:,1),tmp_vec(:,2),tmp_vec(:,3))
