% eLeaf: 3D model of rice leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 1.2.4

function [tri,pts]=ply_read_xy(Path)
%% ply_read_xy
[fid,Msg]=fopen(Path,'rt');
num_pts=fscanf(fid,'%d\n',1);
num_tri=fscanf(fid,'%d\n',1);
pts=zeros(num_pts,3);
tri=zeros(num_tri,3);
for i=1:num_pts
    pts(i,:)=fscanf(fid,'%e %e %e\n',3);
end
for i=1:num_tri
    tri(i,:)=fscanf(fid,'3 %d %d %d\n',3);
end
tri=tri+1;
% trisurf(tri,pts(:,1),pts(:,2),pts(:,3));
% axis equal;
fclose(fid);
end
