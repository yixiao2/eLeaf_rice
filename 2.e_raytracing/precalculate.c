/**************************************************************
*% eLeaf: 3D model of rice leaf photosynthesis
*% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
*% @author: Yi Xiao <yixiao20@outlook.com>
*% @version: 1.2.5
**************************************************************/

/***********************************************
*    version 11/11/11
*
**********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Defs.h"
#include "Global.h"
#include "Trace.h"

void precalculate(Object *obj)
{
	int i;
	Vec v1,v2,v3;
	Triangle *cur_tri;
	Point *p[3];
	
	//start from triangle1, initial 3 edges, initial triangle_info
	int tnum=obj->S->tnum;
	int *tri_info;
	tri_info=(int *)malloc(sizeof(int)*(tnum+1));
	int *edgeunmatch;
	edgeunmatch=(int *)malloc(sizeof(int)*(tnum+1)*6);//3 edges * 2 points per edges
	int *edgeunmatch_info;
	edgeunmatch_info=(int *)malloc(sizeof(int)*(tnum+1)*3);
	
	for(i=0;i<tnum;i++)
	{
		tri_info[i]=0;//0 means this tri is not covered yet
		edgeunmatch_info[i]=0;//0 means this place is not added, 1 means added, 2 means matched
	}
	cur_tri=obj->S->tri+0;	// start from tri0, pt0-1-2
	p[0]=obj->S->vertex+(cur_tri->n[0]);
	p[1]=obj->S->vertex+(cur_tri->n[1]);
	p[2]=obj->S->vertex+(cur_tri->n[2]);	// three vertex of the triangle
	VecSub(p[0]->coordinate,p[1]->coordinate,v1);
	VecSub(p[0]->coordinate,p[2]->coordinate,v2);
	VecCross(v1,v2,v3);
	Normalize(v3);
	VecCopy(v3,cur_tri->Nor);
	
	edgeunmatch_info[0]=1;
	edgeunmatch[0]=cur_tri->n[0];
	edgeunmatch[1]=cur_tri->n[1];
	edgeunmatch_info[1]=1;
	edgeunmatch[2]=cur_tri->n[1];
	edgeunmatch[3]=cur_tri->n[2];
	edgeunmatch_info[2]=1;
	edgeunmatch[4]=cur_tri->n[2];
	edgeunmatch[5]=cur_tri->n[0];//0-1 1-2 2-0
	tri_info[0]=1;//tri0 is done
	
	//search edges, add edge, delete triangle, until all triangle have been covered
	int tri_info_sum,flag_s_edgeunmatch=0,flag_e_edgeunmatch=3,flag_stop=0,tmp_edge,tmp_edge_pt[2],j1,j2,j3,flag_stop2,count_loop=0;
	do
	{
		//search edges from edgeunmatch_info start from No.flag_s_edgeunmatch end edgeunmatch_info=0
		flag_stop=0;
		for(i=flag_s_edgeunmatch;i<flag_e_edgeunmatch&&flag_stop==0;i++)
		{
			if(edgeunmatch_info[i]==1)
			{
				flag_stop=1;
				tmp_edge=i;
				tmp_edge_pt[0]=edgeunmatch[tmp_edge*2];
				tmp_edge_pt[1]=edgeunmatch[tmp_edge*2+1];
			}
		}
		if(flag_stop==0)
		{
			printf("fatal error: disconnected surface\n");
			exit(1);
		}
		//search the edge from tri from tri_info==0
		flag_stop=0;
		for(i=0;i<tnum&&flag_stop==0;i++)
		{
			if(tri_info[i]!=1&&tri_info[i]!=-1)
			{
				cur_tri=obj->S->tri+i;//cur_tri->n[0],n[1],n[2]
				j1=0;
				while(tmp_edge_pt[0]!=cur_tri->n[j1]&&j1<3)
					j1++;
				if(j1!=3)//find a match for one points
				{
					j2=0;
					while(tmp_edge_pt[1]!=cur_tri->n[j2]&&j2<3)
						j2++;
					if(j2!=3)//find a edge
					{
						edgeunmatch_info[tmp_edge]=2;
						//flag_s_edgeunmatch=tmp_edge+1;
						flag_stop=1;
						j3=3-j1-j2;
						//tmp_edge_pt[0]-tmp_edge_pt[1] ==== if j2j1j3 j3j2j1 j1j3j2; then tri_info=1
						if((j2==0&&j1==1)||(j2==1&&j1==2)||(j2==2&&j1==0))
							tri_info[i]=1;
						else
							tri_info[i]=-1;
					}
				}
			}
		}
		if(flag_stop==0)
		{
			printf("fatal error 2\n");
			exit(1);
		}
		//add edges without repetition; based on cur_tri j1 j2 j3;edge j1-j3 j3-j2
		flag_stop=0;//j1-j3
		flag_stop2=0;//j3-j2
		for(i=flag_s_edgeunmatch;i<flag_e_edgeunmatch&&(flag_stop==0||flag_stop2==0);i++)
		{
			if(edgeunmatch_info[i]==1)
			{
				tmp_edge=i;
				tmp_edge_pt[0]=edgeunmatch[tmp_edge*2];
				tmp_edge_pt[1]=edgeunmatch[tmp_edge*2+1];
				
				if(flag_stop==0)
				{
					if(((tmp_edge_pt[0]==cur_tri->n[j1])&&(tmp_edge_pt[1]==cur_tri->n[j3]))||((tmp_edge_pt[0]==cur_tri->n[j3])&&(tmp_edge_pt[1]==cur_tri->n[j1])))
					{
						flag_stop=1;
						edgeunmatch_info[i]=2;
					}
				}
				if(flag_stop2==0)
				{
					if(((tmp_edge_pt[0]==cur_tri->n[j3])&&(tmp_edge_pt[1]==cur_tri->n[j2]))||((tmp_edge_pt[0]==cur_tri->n[j2])&&(tmp_edge_pt[1]==cur_tri->n[j3])))
					{
						flag_stop2=1;
						edgeunmatch_info[i]=2;
					}
				}
			}
		}
		if(flag_stop==0)//add j1-j3 to edgeunmatch
		{
			edgeunmatch[flag_e_edgeunmatch*2]=cur_tri->n[j1];
			edgeunmatch[flag_e_edgeunmatch*2+1]=cur_tri->n[j3];
			edgeunmatch_info[flag_e_edgeunmatch]=1;
			flag_e_edgeunmatch++;
		}
		if(flag_stop2==0)//add j3-j2 to edgeunmatch
		{
			edgeunmatch[flag_e_edgeunmatch*2]=cur_tri->n[j3];
			edgeunmatch[flag_e_edgeunmatch*2+1]=cur_tri->n[j2];
			edgeunmatch_info[flag_e_edgeunmatch]=1;
			flag_e_edgeunmatch++;
		}
		
		tri_info_sum=1;//test stop
		for(i=0;i<tnum;i++)
			tri_info_sum*=tri_info[i];
		count_loop=count_loop+1;
	}while(tri_info_sum==0&&count_loop<tnum);
	
	//volume calculation to determine "in" or "out"
	//mesh centroid
	int pnum=obj->S->pnum;
	Vec center;
	center[0]=0.0;
	center[1]=0.0;
	center[2]=0.0;
	for(i=0;i<pnum;i++)
	{
		center[0]+=(obj->S->vertex+i)->coordinate[0];
		center[1]+=(obj->S->vertex+i)->coordinate[1];
		center[2]+=(obj->S->vertex+i)->coordinate[2];
	}
	center[0]=center[0]/pnum;
	center[1]=center[1]/pnum;
	center[2]=center[2]/pnum;
	//calculate volume_sum
	Vec q[3];
	double volume_sum=0.0,tmp_dbl;
	for(i=0;i<tnum;i++)
	{
		cur_tri=obj->S->tri+i;	// calculate normal vector
		p[0]=obj->S->vertex+(cur_tri->n[0]);
		p[1]=obj->S->vertex+(cur_tri->n[1]);
		p[2]=obj->S->vertex+(cur_tri->n[2]);
		VecSub(center,p[0]->coordinate,q[0]);
		VecSub(center,p[1]->coordinate,q[1]);
		VecSub(center,p[2]->coordinate,q[2]);
		VecCross(q[0],q[1],v1);
		tmp_dbl=VecDot(v1,q[2]);
		if(flag_debug_precal==1)
		{
			if(tmp_dbl*tri_info[i]>0.0)
				printf("+ ");
			else
				printf("- "); 
		}
		volume_sum+=(tmp_dbl)*tri_info[i]/6;
	}
	if(volume_sum<0.0)
	{
		for(i=0;i<tnum;i++)
			tri_info[i]=-tri_info[i];
	}
	
	//finally calculate normal vector
	for(i=0;i<tnum;i++)
	{
		cur_tri=obj->S->tri+i;	// calculate normal vector
		p[0]=obj->S->vertex+(cur_tri->n[0]);
		p[1]=obj->S->vertex+(cur_tri->n[1]);
		p[2]=obj->S->vertex+(cur_tri->n[2]);
		VecSub(p[0]->coordinate,p[1]->coordinate,v1);
		VecSub(p[0]->coordinate,p[2]->coordinate,v2);
		if(tri_info[i]==1)
			VecCross(v1,v2,v3);
		else
			VecCross(v2,v1,v3);
		Normalize(v3);
		VecCopy(v3,cur_tri->Nor);
	}
}

void precalculate_bak(Object *obj)	// obj->S->tri->Nor;
{
	int i;
	Vec v1,v2,v3;
	Triangle *cur_tri;
	Point *p[3];
//	printf("vertexnum %d; trinum %d\n",obj->S->pnum,obj->S->tnum);
	for(i=0;i<obj->S->tnum;i++)
	{
		cur_tri=obj->S->tri+i;	// current triangle
		p[0]=obj->S->vertex+(cur_tri->n[0]);
		p[1]=obj->S->vertex+(cur_tri->n[1]);
		p[2]=obj->S->vertex+(cur_tri->n[2]);	// three vertex of the triangle
		VecSub(p[0]->coordinate,p[1]->coordinate,v1);
		VecSub(p[0]->coordinate,p[2]->coordinate,v2);
		VecCross(v1,v2,v3);
		Normalize(v3);
	//	printf("triangle%d's normal vector:%e %e %e\n",i,v3[0],v3[1],v3[2]);
		
		VecCopy(v3,cur_tri->Nor);
	}
}

