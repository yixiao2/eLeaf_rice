/**************************************************************
*% eLeaf: 3D model of rice leaf photosynthesis
*% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
*% @author: Yi Xiao <yixiao20@outlook.com>
*% @version: 1.2.6
**************************************************************/

/***********************************************
*    version 2017/06/14
*    ray tracing v2 for rice
*    version 2021/10
*    with the same trace.c as eLeaf_dicot
**********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "Defs.h"	// other head files
#include "Global.h"
#include "Trace.h"
#include "precalculate.h"

int main(int argc, char **argv)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * *
	 *		Part1. read in the object and precalculate the normal vector
	 * * * * * * * * * * * * * * * * * * * * * * * * * */
	time_t start=time(NULL);

	SAC_water=atof(*++argv);//Now, input unit m-1, previously, 1e2*atof(*++argv);
	SAC_chlab=atof(*++argv);//Now, input unit m2g-1, previously, 1e-4*atof(*++argv);
	SAC_cyto=SAC_water;
	SAC_air=0.0;
	printf("SAC_water:%e SAC_chlab:%e\n",SAC_water,SAC_chlab);
	
	if(flag_opt_absorbevents==1)
	{
		if((fout_absorbevents=fopen(*++argv,"wb"))==NULL)
		{
			printf("error on open file results_absorbevents_?.txt\n");
		}
	}
	if(flag_RT_file4plot==1)
	{
		if((fout_file4plot=fopen(*++argv,"wb"))==NULL)
		{
			printf("error on open file results_RTfile4plot_?.txt\n");
		}
	}
	
	FILE *fp;
	char cell[25];
	char chl[22];
	char vacuole[23];
	char str_num[3];

	int i,j,loop1;
	int dim;
	
	Object *curobj;
	Point *temp_p=NULL;
	Triangle *temp_t=NULL;	
	int vertex_num,tri_num;
	int tmp_bsc_num,tmp_msc_num;
	//char tmp_buf[128];
	
	//////////////////////read in count_chl///////////////
	fp=fopen(*++argv,"r");
	fscanf(fp,"%d %d\n",&tmp_bsc_num,&tmp_msc_num);
	for(loop1=1;loop1<=(msall_num);loop1++)
	{
		fscanf(fp,"%d %lf\n",&(ms_chl_num[loop1-1]),&(ms_chl_con[loop1-1]));
	}
	fclose(fp);
	
	//////////////////////////read in MS///////////////////////	
	printf("read in MS\n");
	count_ms=0;	
	for(loop1=1;loop1<=(msall_num);loop1++)
	{	   
		count_ms++;
		/********************
		*    cell
		********************/
		strcpy(cell,"MS/ms");
		sprintf(str_num,"%d",loop1);
		strcat(cell,str_num);
		strcat(cell,".ply");
		
		if((fp=fopen(cell,"r"))==NULL)
		{
			printf("error on open cell: %s\n",cell);
		}
		fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//1st line -- ply
		//fgets(tmp_buf, sizeof tmp_buf, fp);//2nd line -- format ascii 1.0
		//fgets(tmp_buf, sizeof tmp_buf, fp);//3rd line -- element vertex 29
		//sscanf(tmp_buf,"%*s %*s %d", &vertex_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//4th line -- property float x
		//fgets(tmp_buf, sizeof tmp_buf, fp);//5th line -- property float y
		//fgets(tmp_buf, sizeof tmp_buf, fp);//6th line -- property float z
		//fgets(tmp_buf, sizeof tmp_buf, fp);//7th line -- element face 54
		//sscanf(tmp_buf,"%*s %*s %d", &tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//8th line -- property list uchar int vertex_indices
		//fgets(tmp_buf, sizeof tmp_buf, fp);//9th line -- end_header

		curobj=(Object *)malloc(sizeof(Object));
		curobj->S=(Surface *)malloc(sizeof(Surface));
		curobj->S->pnum=vertex_num;
		curobj->S->tnum=tri_num;
		curobj->S->vertex=(Point *)malloc(sizeof(Point)*vertex_num);
		curobj->S->tri=(Triangle *)malloc(sizeof(Triangle)*tri_num);
		temp_p=curobj->S->vertex;
		temp_t=curobj->S->tri;
		for(i=0;i<vertex_num;i++)
		{
			fscanf(fp,"%lf %lf %lf\n",&temp_p->coordinate[0],&temp_p->coordinate[1],&temp_p->coordinate[2]);		
			temp_p+=1;
		}
		for(i=0;i<tri_num;i++)
		{
			fscanf(fp,"%d %d %d %d\n",&dim,&temp_t->n[0],&temp_t->n[1],&temp_t->n[2]);
			temp_t+=1;
		}
		fclose(fp);
		precalculate(curobj);
		p_cell_ms[count_ms-1]=curobj;
				
		/********************
		*    chl
		********************/
		for(i=1;i<=ms_chl_num[count_ms-1];i++)
		{
			strcpy(chl,"MS/ms");
			sprintf(str_num,"%d",loop1);
			strcat(chl,str_num);
			sprintf(str_num,"%d",i);
			strcat(chl,"c");
			strcat(chl,str_num);
			strcat(chl,".ply");
			//printf("chl:%s\n",chl);

			if((fp=fopen(chl,"r"))==NULL)
			{
				printf("error on open %s\n",chl);
			}

			fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
			//fgets(tmp_buf, sizeof tmp_buf, fp);//1st line -- ply
			//fgets(tmp_buf, sizeof tmp_buf, fp);//2nd line -- format ascii 1.0
			//fgets(tmp_buf, sizeof tmp_buf, fp);//3rd line -- element vertex 29
			//sscanf(tmp_buf,"%*s %*s %d", &vertex_num);
			//fgets(tmp_buf, sizeof tmp_buf, fp);//4th line -- property float x
			//fgets(tmp_buf, sizeof tmp_buf, fp);//5th line -- property float y
			//fgets(tmp_buf, sizeof tmp_buf, fp);//6th line -- property float z
			//fgets(tmp_buf, sizeof tmp_buf, fp);//7th line -- element face 54
			//sscanf(tmp_buf,"%*s %*s %d", &tri_num);
			//fgets(tmp_buf, sizeof tmp_buf, fp);//8th line -- property list uchar int vertex_indices
			//fgets(tmp_buf, sizeof tmp_buf, fp);//9th line -- end_header

			curobj=(Object *)malloc(sizeof(Object));
			curobj->S=(Surface *)malloc(sizeof(Surface));
			curobj->S->pnum=vertex_num;
			curobj->S->tnum=tri_num;
			curobj->S->vertex=(Point *)malloc(sizeof(Point)*vertex_num);
			curobj->S->tri=(Triangle *)malloc(sizeof(Triangle)*tri_num);
			temp_p=curobj->S->vertex;
			temp_t=curobj->S->tri;
			for(j=0;j<vertex_num;j++)
			{
				fscanf(fp,"%lf %lf %lf\n",&temp_p->coordinate[0],&temp_p->coordinate[1],&temp_p->coordinate[2]);
				temp_p+=1;
			}
			for(j=0;j<tri_num;j++)
			{
				fscanf(fp,"%d %d %d %d\n",&dim,&temp_t->n[0],&temp_t->n[1],&temp_t->n[2]);
				temp_t+=1;
			}
			fclose(fp);
			precalculate(curobj);
			p_chl_ms[count_ms-1][i-1]=curobj;
		}
		
		/***********************
		*    vacuole
		***********************/
		strcpy(vacuole,"MS/ms");
		sprintf(str_num,"%d",loop1);
		strcat(vacuole,str_num);
		strcat(vacuole,"v.ply");
		
		if((fp=fopen(vacuole,"r"))==NULL)
		{
			printf("error on open vacuole: %s\n",vacuole);
		}
		fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//1st line -- ply
		//fgets(tmp_buf, sizeof tmp_buf, fp);//2nd line -- format ascii 1.0
		//fgets(tmp_buf, sizeof tmp_buf, fp);//3rd line -- element vertex 29
		//sscanf(tmp_buf,"%*s %*s %d", &vertex_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//4th line -- property float x
		//fgets(tmp_buf, sizeof tmp_buf, fp);//5th line -- property float y
		//fgets(tmp_buf, sizeof tmp_buf, fp);//6th line -- property float z
		//fgets(tmp_buf, sizeof tmp_buf, fp);//7th line -- element face 54
		//sscanf(tmp_buf,"%*s %*s %d", &tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//8th line -- property list uchar int vertex_indices
		//fgets(tmp_buf, sizeof tmp_buf, fp);//9th line -- end_header

		curobj=(Object *)malloc(sizeof(Object));
		curobj->S=(Surface *)malloc(sizeof(Surface));
		curobj->S->pnum=vertex_num;
		curobj->S->tnum=tri_num;
		curobj->S->vertex=(Point *)malloc(sizeof(Point)*vertex_num);
		curobj->S->tri=(Triangle *)malloc(sizeof(Triangle)*tri_num);
		temp_p=curobj->S->vertex;
		temp_t=curobj->S->tri;
		for(i=0;i<vertex_num;i++)
		{
			fscanf(fp,"%lf %lf %lf\n",&temp_p->coordinate[0],&temp_p->coordinate[1],&temp_p->coordinate[2]);		
			temp_p+=1;
		}
		for(i=0;i<tri_num;i++)
		{
			fscanf(fp,"%d %d %d %d\n",&dim,&temp_t->n[0],&temp_t->n[1],&temp_t->n[2]);
			temp_t+=1;
		}
		fclose(fp);
		precalculate(curobj);
		p_vac_ms[count_ms-1]=curobj;
	}
	printf("\n");

	//////////////////////////read in leaf///////////////////////
	printf("read in the leaf\n");
	if((fp=fopen("leaf","r"))==NULL)
	{
		printf("error on open file leaf\n");
	}
	fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
	//fgets(tmp_buf, sizeof tmp_buf, fp);//1st line -- ply
	//fgets(tmp_buf, sizeof tmp_buf, fp);//2nd line -- format ascii 1.0
	//fgets(tmp_buf, sizeof tmp_buf, fp);//3rd line -- element vertex 29
	//sscanf(tmp_buf,"%*s %*s %d", &vertex_num);
	//fgets(tmp_buf, sizeof tmp_buf, fp);//4th line -- property float x
	//fgets(tmp_buf, sizeof tmp_buf, fp);//5th line -- property float y
	//fgets(tmp_buf, sizeof tmp_buf, fp);//6th line -- property float z
	//fgets(tmp_buf, sizeof tmp_buf, fp);//7th line -- element face 54
	//sscanf(tmp_buf,"%*s %*s %d", &tri_num);
	//fgets(tmp_buf, sizeof tmp_buf, fp);//8th line -- property list uchar int vertex_indices
	//fgets(tmp_buf, sizeof tmp_buf, fp);//9th line -- end_header

	curobj=(Object *)malloc(sizeof(Object));
	curobj->S=(Surface *)malloc(sizeof(Surface));
	curobj->S->pnum=vertex_num;
	curobj->S->tnum=tri_num;
	curobj->S->vertex=(Point *)malloc(sizeof(Point)*vertex_num);
	curobj->S->tri=(Triangle *)malloc(sizeof(Triangle)*tri_num);
	temp_p=curobj->S->vertex;
	temp_t=curobj->S->tri;
	for(j=0;j<vertex_num;j++)
	{
		fscanf(fp,"%lf %lf %lf\n",&temp_p->coordinate[0],&temp_p->coordinate[1],&temp_p->coordinate[2]);
		temp_p+=1;
	}
	
	for(j=0;j<tri_num;j++)
	{
		fscanf(fp,"%d %d %d %d\n",&dim,&temp_t->n[0],&temp_t->n[1],&temp_t->n[2]);
		temp_t+=1;
	}
	fclose(fp);
	precalculate(curobj);
	p_leaf=curobj;
	
	////////////////////////////read non MS///////////////////////
	printf("read in nonMS cell\n");
	count_nonms=0;	
	//EPL_l
	for(loop1=1;loop1<=nonms_num;loop1++)//upper and lower eppidermis, bulliform cells
	{
		count_nonms++;
		//////////read in cell/////////
		strcpy(cell,"nonMS/ns");
		sprintf(str_num,"%d",loop1);
		strcat(cell,str_num);		
		strcat(cell,".ply");
		//printf("cell:%s\n",cell);

		if((fp=fopen(cell,"r"))==NULL)
		{
			printf("error on open %s\n",cell);
		}		
		fscanf(fp,"%d\n%d\n",&vertex_num,&tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//1st line -- ply
		//fgets(tmp_buf, sizeof tmp_buf, fp);//2nd line -- format ascii 1.0
		//fgets(tmp_buf, sizeof tmp_buf, fp);//3rd line -- element vertex 29
		//sscanf(tmp_buf,"%*s %*s %d", &vertex_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//4th line -- property float x
		//fgets(tmp_buf, sizeof tmp_buf, fp);//5th line -- property float y
		//fgets(tmp_buf, sizeof tmp_buf, fp);//6th line -- property float z
		//fgets(tmp_buf, sizeof tmp_buf, fp);//7th line -- element face 54
		//sscanf(tmp_buf,"%*s %*s %d", &tri_num);
		//fgets(tmp_buf, sizeof tmp_buf, fp);//8th line -- property list uchar int vertex_indices
		//fgets(tmp_buf, sizeof tmp_buf, fp);//9th line -- end_header

		curobj=(Object *)malloc(sizeof(Object));
		curobj->S=(Surface *)malloc(sizeof(Surface));
		curobj->S->pnum=vertex_num;
		curobj->S->tnum=tri_num;
		curobj->S->vertex=(Point *)malloc(sizeof(Point)*vertex_num);
		curobj->S->tri=(Triangle *)malloc(sizeof(Triangle)*tri_num);
		temp_p=curobj->S->vertex;
		temp_t=curobj->S->tri;
		for(i=0;i<vertex_num;i++)
		{
			fscanf(fp,"%lf %lf %lf\n",&temp_p->coordinate[0],&temp_p->coordinate[1],&temp_p->coordinate[2]);		
			temp_p+=1;
		}
		for(i=0;i<tri_num;i++)
		{
			fscanf(fp,"%d %d %d %d\n",&dim,&temp_t->n[0],&temp_t->n[1],&temp_t->n[2]);
			temp_t+=1;
		}
		fclose(fp);
		precalculate(curobj);
		p_cell_ns[count_nonms-1]=curobj;
		printf("finish read in %s\n",cell);
	}
	printf("\n");
		
	printf("finish read in \n");

	/*///////for debug the new precalculate()//////////////
	curobj=p_chl_ms[1][0];//p_leaf p_cell_ms[0] p_chl_ms[0][0] p_vac_ms[0] p_cell_ns[0]
	tri_num=curobj->S->tnum;
	if(flag_debug_precal==1)
	{
		FILE *f_tmp;
		Triangle *tri_tmp;
		f_tmp=fopen("debug_norvec.txt","wb");
		for(i=0;i<tri_num;i++)
		{
			tri_tmp=curobj->S->tri+i;
			fprintf(f_tmp,"%e %e %e\n",tri_tmp->Nor[0],tri_tmp->Nor[1],tri_tmp->Nor[2]);
		}
		fclose(f_tmp);
	}
	///////output normal vector of object///////////////*/
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *			Part2. initialize the object
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	///////////////////////leaf///////////////// 
	curobj=p_leaf;
	curobj->refr_index_o=0.0001;//use refr_index_o here to mirror rays back. Total reflection.
	curobj->refr_index_i=air_refr_index;
	curobj->refr_index_s=-1.0;//-1 means inactive
	curobj->sac_o=0.0;
	curobj->sac_i=0.0;
	curobj->I_absorb_o=0.0;
	curobj->I_absorb_i=0.0;
	for(j=0;j<curobj->S->tnum;j++)
	{
		temp_t=curobj->S->tri+j;
		temp_t->I_absorb_tri_o=0.0;
		temp_t->I_absorb_tri_i=0.0;
		temp_t->count_hit_tri=0;
	}
	curobj->subobj=p_cell_ms[0];
	curobj->nextobj=NULL;
	curobj->belongobj=NULL;
	curobj->celltype=0;
	curobj->cellname=0;
	curobj->chlname=0;
	curobj->vacuolename=0;
	curobj->count_hit=0;
	///////////////////////ms///////////////// 
	for(loop1=1;loop1<=msall_num;loop1++)
	{
		//////// p_cell_ms/////////
		curobj=p_cell_ms[loop1-1];
		curobj->refr_index_o=air_refr_index;
		curobj->refr_index_i=cyto_refr_index;
		curobj->refr_index_s=wall_refr_index;//cell wall is considered specifically in trace.c
		curobj->sac_o=SAC_air;
		curobj->sac_i=SAC_cyto;
		curobj->I_absorb_o=0.0;
		curobj->I_absorb_i=0.0;
		for(j=0;j<curobj->S->tnum;j++)
		{
			temp_t=curobj->S->tri+j;
			temp_t->I_absorb_tri_o=0.0;
			temp_t->I_absorb_tri_i=0.0;
			temp_t->count_hit_tri=0;
		}
		curobj->subobj=p_chl_ms[loop1-1][0];
		if(loop1!=msall_num)
		{
			curobj->nextobj=p_cell_ms[loop1];
		}
		else
		{
			curobj->nextobj=p_cell_ns[0];
		}
		curobj->belongobj=p_leaf;
		curobj->celltype=1;
		curobj->cellname=loop1;
		curobj->chlname=0;
		curobj->vacuolename=0;
		curobj->count_hit=0;
		
		////////////////p_chl_ms////////////////
		for(i=1;i<=ms_chl_num[loop1-1];i++)
		{
			curobj=p_chl_ms[loop1-1][i-1];
			curobj->refr_index_o=cyto_refr_index; // refraction index
			curobj->refr_index_i=chlo_refr_index;
			curobj->refr_index_s=-1.0;// inactive
			
			curobj->sac_o=SAC_cyto;
			curobj->sac_i=SAC_chlab*ms_chl_con[loop1-1]; // same concentration for all chloroplasts
			curobj->I_absorb_o=0.0;
			curobj->I_absorb_i=0.0;
			for(j=0;j<curobj->S->tnum;j++)
			{
				temp_t=curobj->S->tri+j;
				temp_t->I_absorb_tri_o=0.0;
				temp_t->I_absorb_tri_i=0.0;
				temp_t->count_hit_tri=0;
			}
			curobj->subobj=NULL;
			if(i==ms_chl_num[loop1-1])
				curobj->nextobj=p_vac_ms[loop1-1];
			else
				curobj->nextobj=p_chl_ms[loop1-1][i];
			curobj->belongobj=p_cell_ms[loop1-1];
			curobj->celltype=2;
			curobj->cellname=loop1;
			curobj->chlname=i;
			curobj->vacuolename=0;
			curobj->count_hit=0;
		}
		////////////////p_vac_ms////////////////
		curobj=p_vac_ms[loop1-1];
		curobj->refr_index_o=cyto_refr_index; // refraction index
		curobj->refr_index_i=vacu_refr_index;
		curobj->refr_index_s=-1.0;// inactive
		
		curobj->sac_o=SAC_cyto;
		curobj->sac_i=SAC_water;
		curobj->I_absorb_o=0.0;
		curobj->I_absorb_i=0.0;
		for(j=0;j<curobj->S->tnum;j++)
		{
			temp_t=curobj->S->tri+j;
			temp_t->I_absorb_tri_o=0.0;
			temp_t->I_absorb_tri_i=0.0;
			temp_t->count_hit_tri=0;
		}
		curobj->subobj=NULL;
		curobj->nextobj=NULL;
		curobj->belongobj=p_cell_ms[loop1-1];
		curobj->celltype=4;
		curobj->cellname=loop1;
		curobj->chlname=0;
		curobj->vacuolename=1;
		curobj->count_hit=0;
	}
	printf("finish initial the ms\n");

	//////////////////nonMS///////////////////
	for(loop1=1;loop1<=nonms_num;loop1++)
	{
		curobj=p_cell_ns[loop1-1];
		curobj->refr_index_o=air_refr_index;
		curobj->refr_index_i=cyto_refr_index;
		curobj->refr_index_s=wall_refr_index;
		
		curobj->sac_o=SAC_cyto;
		curobj->sac_i=SAC_cyto;
		curobj->I_absorb_o=0.0;
		curobj->I_absorb_i=0.0;
		for(j=0;j<curobj->S->tnum;j++)
		{
			temp_t=curobj->S->tri+j;
			temp_t->I_absorb_tri_o=0.0;
			temp_t->I_absorb_tri_i=0.0;
			temp_t->count_hit_tri=0;
		}
		curobj->subobj=NULL;

		if(loop1!=nonms_num)
		{
			curobj->nextobj=p_cell_ns[loop1];
		}
		else
		{
			curobj->nextobj=NULL;
		}
		
		curobj->belongobj=p_leaf;
		curobj->celltype=3;
		curobj->cellname=loop1;
		curobj->chlname=0;
		curobj->vacuolename=0;
		curobj->count_hit=0;
	}	
	printf("finish initialization\n");
	
	I_discard=0.0;
	I_discard_Rf=0.0;
	I_discard_Tr=0.0;

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *			Part3. generate light and trace
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	//int cutnum_x=500,cutnum_y=25;
	int cutnum_x,cutnum_y;
	int ray_i_start,ray_i_end;
	int ray_j_start,ray_j_end;

	cutnum_x=atoi(*++argv);
	cutnum_y=atoi(*++argv);
	ray_i_start=atoi(*++argv);
	ray_i_end=atoi(*++argv);
	ray_j_start=atoi(*++argv);
	ray_j_end=atoi(*++argv);
	double dx=(xmax-xmin)/(cutnum_x+1),dy=(zmax-zmin)/(cutnum_y+1);
/*		
	Ray *curray;
	curray=(Ray *)malloc(sizeof(Ray));
	curray->P[0]=-2.400000e-05;
	curray->P[1]=-2.313500e-06;
	curray->P[2]=zmax;
	curray->D[0]=-4.114551e-01;
	curray->D[1]=8.075262e-01;
	curray->D[2]=-4.226183e-01;
	curray->I=1.0;
//	curray->interflag=0;
	curray->level=1;
	curray->P_obj=p_leaf;
	curray->interflag=0;
//	curray->P_tri=175;
	printf("\n");
	printf("start ray tracing\n");
	printf("\n");
	
	trace(curray,p_cell_s_o[0]);
*/

	if(flag_randseed==1)
	{
		srand(time(NULL));
	}

	Ray *curray;
	double alpha,beta;

	for(ray_i=ray_i_start;ray_i<ray_i_end;ray_i++) ///// 0 < cutnum
	{
		for(ray_j=ray_j_start;ray_j<ray_j_end;ray_j++)
		{
			curray=(Ray *)malloc(sizeof(Ray));
			// initialize the current light
			debugI=0.0;
			num_chl_hit=0;
			
			curray->P[0]=xmin+(ray_i+1)*dx;
			curray->P[1]=ymax-0.05e-6;
			curray->P[2]=zmin+(ray_j+1)*dy;
		
			curray->D[0]=0.0;
			curray->D[1]=-1.0;
			curray->D[2]=0.0;
			curray->P_tri=4;//previously -1; but cause issue to trace_recal() based on absorbevent files
						
			if(flag_pertube_dir==1)
			{
				//alpha=rand()%178+1;
				alpha=(double)rand()/RAND_MAX;//0-1 degree
				//printf("alpha:%e\n",alpha);
				curray->D[1]=-cos(alpha/180*PI);
				beta=rand()/360;
				curray->D[0]=sin(alpha/180*PI)*sin(beta/180*PI);
				curray->D[2]=sin(alpha/180*PI)*cos(beta/180*PI);
			}

			curray->I=1.0; 
			curray->interflag=0;
			curray->level=1;
			curray->P_obj=p_leaf;
			curray->trace_obj=p_leaf->subobj;
			curray->refr_index=air_refr_index;
			curray->outin_P=1;
			
			printf("\n************BEGIN************\n");
			printf("*** Trace ray i=%d j=%d ***\n",ray_i,ray_j);
			//printf("P:%e %e %e\n",curray->P[0],curray->P[1],curray->P[2]);
			//printf("D:%e %e %e\n",curray->D[0],curray->D[1],curray->D[2]);
			if(flag_RT_file4plot==1)
			{
				fprintf(fout_file4plot,"%d %d %e %e %e\n",ray_i,ray_j,curray->P[0],curray->P[1],curray->P[2]);
			}
			trace(curray,p_leaf->subobj);
			
			sum_p();
		}
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *			Part4. output file
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	int k;
	FILE *fout_tri, *fout_srf;
	double sum_I_absorb_all=0.0;
	
	if(flag_output_results_tri!=0)
	{
		if((fout_tri=fopen(*++argv,"wb"))==NULL)
		{
			printf("error on open file result tri\n");
		}
	}
	if(flag_output_results_srf!=0)
	{
		if((fout_srf=fopen(*++argv,"wb"))==NULL)
		{
			printf("error on open file result cell\n");
		}
	}

	//////////// MS ///////////////////
	for(i=1;i<=msall_num;i++)
	{
		curobj=p_cell_ms[i-1];
		sum_I_absorb_all+=(curobj->I_absorb_o+curobj->I_absorb_i);
		if(flag_output_results_tri!=0)
		{
			for(k=0;k<curobj->S->tnum;k++)
			{
				fprintf(fout_tri,"%d %d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,(curobj->S->tri+k)->count_hit_tri,
										(curobj->S->tri+k)->I_absorb_tri_o,(curobj->S->tri+k)->I_absorb_tri_i);
			}
		}
		if(flag_output_results_srf!=0)
		{
			fprintf(fout_srf,"%d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb_o,curobj->I_absorb_i);
		}
		
		for(j=1;j<=ms_chl_num[i-1];j++)
		{
			curobj=p_chl_ms[i-1][j-1];
			sum_I_absorb_all+=(curobj->I_absorb_o+curobj->I_absorb_i);
			if(flag_output_results_tri!=0)
			{
				for(k=0;k<curobj->S->tnum;k++)
				{
					fprintf(fout_tri,"%d %d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,
											(curobj->S->tri+k)->count_hit_tri,(curobj->S->tri+k)->I_absorb_tri_o,(curobj->S->tri+k)->I_absorb_tri_i);
				}
			}
			if(flag_output_results_srf!=0)
			{
				fprintf(fout_srf,"%d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb_o,curobj->I_absorb_i);
			}
		}
		
		curobj=p_vac_ms[i-1];
		sum_I_absorb_all+=(curobj->I_absorb_o+curobj->I_absorb_i);
		if(flag_output_results_tri!=0)
		{
			for(k=0;k<curobj->S->tnum;k++)
			{
				fprintf(fout_tri,"%d %d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,
										(curobj->S->tri+k)->count_hit_tri,(curobj->S->tri+k)->I_absorb_tri_o,(curobj->S->tri+k)->I_absorb_tri_i);
			}
		}
		if(flag_output_results_srf!=0)
		{
			fprintf(fout_srf,"%d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb_o,curobj->I_absorb_i);
		}
	}
	
	//////////// nonMS ///////////////////
	for(i=1;i<=nonms_num;i++)
	{
		curobj=p_cell_ns[i-1];
		sum_I_absorb_all+=(curobj->I_absorb_o+curobj->I_absorb_i);
		if(flag_output_results_tri!=0)
		{
			for(k=0;k<curobj->S->tnum;k++)
			{
				fprintf(fout_tri,"%d %d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,
										(curobj->S->tri+k)->count_hit_tri,(curobj->S->tri+k)->I_absorb_tri_o,(curobj->S->tri+k)->I_absorb_tri_i);
			}
		}
		if(flag_output_results_srf!=0)
		{
			fprintf(fout_srf,"%d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb_o,curobj->I_absorb_i);	
		}
	}
	if(flag_output_results_tri!=0)
		fclose(fout_tri);
	if(flag_output_results_srf!=0)
		fclose(fout_srf);
	
	FILE *fout_sum;
	if((fout_sum=fopen(*++argv,"wb"))==NULL)
	{
		printf("error open file result_sum\n");
	}
	fprintf(fout_sum,"Absorption %e\n",sum_I_absorb_all);
	fprintf(fout_sum,"Reflectance %e\n",I_discard_Rf);
	fprintf(fout_sum,"Transmitance %e\n",I_discard_Tr);
	fprintf(fout_sum,"Discard %e\n",I_discard);
	fprintf(fout_sum,"Sum_I %e\n",sum_I_absorb_all+I_discard_Rf+I_discard_Tr+I_discard);
	
	time_t end=time(NULL);
	printf("time:%f sec\n",difftime(end,start));
	fprintf(fout_sum,"Time %f\n",difftime(end,start));
	fclose(fout_sum);
	
	if(flag_opt_absorbevents==1)
	{
		fclose(fout_absorbevents);
	}
	if(flag_RT_file4plot==1)
	{
		fclose(fout_file4plot);
	}

	return 1;
}
