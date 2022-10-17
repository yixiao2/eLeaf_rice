/**************************************************************
*% eLeaf: 3D model of rice leaf photosynthesis
*% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
*% @author: Yi Xiao <yixiao20@outlook.com>
*% @version: 1.2.5
**************************************************************/

/***********************************************
*    version 2021 July
*    read results_absorbevents_* and recalculate RT_path.
*    output a ab_tri file, a ab_sum file the absorption profile
*
*    modified from inner_check.c
*
*    Command:
*    ./trace_recal new_SAC_water new_SAC_chl count_chl4RT_new cutnum_x cutnum_y 
*                  cutnum_x_perfile cutnum_y_perfile
*                  absorbevent_filename_prefix absorbsum_filename_prefix
*                  num_layers
*                  fout_abtri fout_absrf fout_abprofile fout_rtsum
*
*    Demo:
*    ./trace_recal 0.0075 0.0610 count_chl4RT 500 500 
*                  5 5 
*                  results_absorbevents_450nm_500x_ results_sum_450nm_500x_ 
*                  10 
*                  test_abtri test_absrf test_abprofile test_rtsum
*
*    trace_recal can also be used to merge results file if use the SAME 
*                  - SAC_water, 
*                  - SAC_chl, 
*                  - count_chl4RT
**********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "Defs.h"	// other head files
#include "Global.h"
//#include "Trace.h"
#include "precalculate.h"

int main(int argc, char **argv)
{
	/* * * * * * * * * * * * * * * * * * * * * * * * * *
	 *		Part1. read in the object and initialize
	 * * * * * * * * * * * * * * * * * * * * * * * * * */
	time_t start=time(NULL);
	
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

	SAC_water=atof(*++argv);//Now, input unit m-1, previously, 1e2*atof(*++argv);
	SAC_chlab=atof(*++argv);//Now, input unit m2g-1, previously, 1e-4*atof(*++argv);
	SAC_cyto=SAC_water;
	SAC_air=0.0;
	
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
		//printf("finish read in %s\n",cell);
	}
	printf("\n");
		
	printf("finish read in all ply files\n");
	
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 *	Part2. initialize the object
	 *             reinitialize with new SAC and [chl]
	 *             new [chl] calculated from new partition ratio with a MATLAB script
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	///////////////////////leaf///////////////// 
	curobj=p_leaf;
	//curobj->refr_index_o=0.0001;//use refr_index_o here to mirror rays back. Total reflection.
	//curobj->refr_index_i=air_refr_index;
	//curobj->refr_index_s=-1.0;//-1 means inactive
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
		//curobj->refr_index_o=air_refr_index;
		//curobj->refr_index_i=cyto_refr_index;
		//curobj->refr_index_s=wall_refr_index;//cell wall is considered specifically in trace.c
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
			//curobj->refr_index_o=cyto_refr_index; // refraction index
			//curobj->refr_index_i=chlo_refr_index;
			//curobj->refr_index_s=-1.0;// inactive
			
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
		//curobj->refr_index_o=cyto_refr_index; // refraction index
		//curobj->refr_index_i=vacu_refr_index;
		//curobj->refr_index_s=-1.0;// inactive
		
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
		//curobj->refr_index_o=air_refr_index;
		//curobj->refr_index_i=cyto_refr_index;
		//curobj->refr_index_s=wall_refr_index;
		
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

	/*****************************************************************
	*       load results_absorbevents_* and recalculate RT_path
	*****************************************************************/
	Object *tmpobj;	
	int loop_file,count_fail=0;
	int ray_i,ray_j,P_celltype,P_cellname,P_chlname,P_vacuolename,P_tri,P_outin_P,ray_status;
	int ray_i0,ray_j0,flag_recal;
	double dis,tmp_I,cur_I;
	char filename_abevent[100],tmp_filename_abevent[100],filename_absum[100],tmp_filename_absum[100], tmp_str[10];

	int cutnum_x,cutnum_y;
	cutnum_x=atoi(*++argv);
	cutnum_y=atoi(*++argv);
	int cutnum_x_perfile,cutnum_y_perfile;
	cutnum_x_perfile=atoi(*++argv);
	cutnum_y_perfile=atoi(*++argv);
	int max_file_idx;
	max_file_idx=(cutnum_x/cutnum_x_perfile)*(cutnum_y/cutnum_y_perfile);

	I_discard=0.0;//defined in defs.h
	I_discard_Rf=0.0;
	I_discard_Tr=0.0;

	strcpy(filename_abevent,*++argv);//readin the suffix of absorbevent files, e.g "results_absorbevents_450nm_500x_"
	strcpy(filename_absum,*++argv);//readin the suffix of absorbevent files, e.g "results_sum_450nm_500x_"

	for(loop_file=1;loop_file<=max_file_idx;loop_file++)
	{
		//strcpy(file_name,"results_absorbevents_450nm_500x_");
		sprintf(tmp_filename_abevent,"%s",filename_abevent);
		sprintf(tmp_filename_absum,"%s",filename_absum);
		sprintf(tmp_str,"%d",loop_file);
		strcat(tmp_filename_abevent, tmp_str);
		strcat(tmp_filename_absum, tmp_str);

		if((fp=fopen(tmp_filename_absum,"r"))==NULL)
		{
			// printf("Sum File does not exist: %s\n",tmp_filename_absum);
			count_fail=count_fail+1;
		}
		else
		{
			fclose(fp);//close fp for sum files
			if((fp=fopen(tmp_filename_abevent,"r"))==NULL)
			{
				printf("[Error] Sum File exists, but Abevent File does not: %s\n",tmp_filename_abevent);
				exit(1);
			}

			ray_i0=-1;//(int)(floor(loop_file/(cutnum_y/cutnum_y_perfile)))*cutnum_x_perfile;
			ray_j0=-1;//(loop_file%(cutnum_x/cutnum_x_perfile))*cutnum_y_perfile-cutnum_y_perfile;
			cur_I=0.0;
			flag_recal=0;	

			// additional %d at last represents ray status, 0=absorb, 1=reflect, 2=transmit
			while(fscanf(fp,"%d %d %lf %d %d %d %d %d %d %lf %d\n",&ray_i,&ray_j,&dis,&P_celltype,&P_cellname,&P_chlname,&P_vacuolename,&P_tri,&P_outin_P,&tmp_I,&ray_status)!=EOF)
			{
				if(ray_i!=ray_i0||ray_j!=ray_j0)
				{
					// check cur_I which is actually the light intensity of last ray
					if(cur_I>discardI)
					{
						printf("ErrorMsg: discard light intensity > discardI\n");
						printf("    ray_i=%d, ray_j=%d.\n",ray_i0,ray_j0);
						printf("    ray_I=%e\n",cur_I);
						exit(1);
					}

					ray_i0=ray_i;
					ray_j0=ray_j;
					cur_I=tmp_I;//should = 1.0;
					flag_recal=1;
					if(ray_status==1)
					{
						I_discard_Rf=I_discard_Rf+cur_I;
						cur_I=0.0;
					}
					else if(ray_status==2)
					{
						I_discard_Tr=I_discard_Tr+cur_I;
						cur_I=0.0;
					}
					else//ray_status==0
					//inactive if initial ray is from leaf_box, but possible for future usage
					{
						//[P_celltype,P_cellname,P_chlname,P_vacuolename]->ray->P_obj?
						//celltype: 0=leaf, 1=ms cellwall, 2=ms chl, 3=nonms cellwall, 4=ms vac
						switch(P_celltype)
						{
							case 0:
								tmpobj=p_leaf;
								break;
							case 1:
								tmpobj=p_cell_ms[P_cellname-1];
								break;
							case 2:
								tmpobj=p_chl_ms[P_cellname-1][P_chlname-1];
								break;
							case 3:
								tmpobj=p_cell_ns[P_cellname-1];
								break;
							case 4:
								tmpobj=p_vac_ms[P_cellname-1];
								break;									
						}
						if(P_outin_P==1)
						{
							////P of ray is on the inner side of surface, i.e. ray is from the mother node
							tmp_I=cur_I*exp(-(tmpobj->sac_i)*(dis));
							(tmpobj->S->tri+P_tri)->I_absorb_tri_i+=(cur_I-tmp_I);
							tmpobj->I_absorb_i+=(cur_I-tmp_I);
						}
						else
						{
							////P of ray is on the outter side of surface, i.e. ray is reflected from a neighbour node
							tmp_I=cur_I*exp(-(tmpobj->sac_o)*(dis));
							(tmpobj->S->tri+P_tri)->I_absorb_tri_o+=(cur_I-tmp_I);
							tmpobj->I_absorb_o+=(cur_I-tmp_I);
						}
						cur_I=tmp_I;
					}
					if(cur_I<discardI)
					{
						flag_recal=0;
						I_discard=I_discard+cur_I;
					}
				}
				else
				{
					if(flag_recal==1)//flag_recal becomes 0 if cur_I < discard_I for a specific ray
					{
						if(ray_status==1)
						{
							I_discard_Rf=I_discard_Rf+cur_I;
							cur_I=0.0;
						}
						else if(ray_status==2)
						{
							I_discard_Tr=I_discard_Tr+cur_I;
							cur_I=0.0;
						}
						else//ray_status==0
						{
							//[P_celltype,P_cellname,P_chlname,P_vacuolename]->ray->P_obj?
							//celltype: 0=leaf, 1=ms cellwall, 2=ms chl, 3=nonms cellwall, 4=ms vac
							switch(P_celltype)
							{
								case 0:
									tmpobj=p_leaf;
									break;
								case 1:
									tmpobj=p_cell_ms[P_cellname-1];
									break;
								case 2:
									tmpobj=p_chl_ms[P_cellname-1][P_chlname-1];
									break;
								case 3:
									tmpobj=p_cell_ns[P_cellname-1];
									break;
								case 4:
									tmpobj=p_vac_ms[P_cellname-1];
									break;									
							}
							

							if(P_outin_P==1)
							{
								////P of ray is on the inner side of surface, i.e. ray is from the mother node
								tmp_I=cur_I*exp(-(tmpobj->sac_i)*(dis));
								(tmpobj->S->tri+P_tri)->I_absorb_tri_i+=(cur_I-tmp_I);
								tmpobj->I_absorb_i+=(cur_I-tmp_I);
							}
							else
							{
								////P of ray is on the outter side of surface, i.e. ray is reflected from a neighbour node
								tmp_I=cur_I*exp(-(tmpobj->sac_o)*(dis));
								(tmpobj->S->tri+P_tri)->I_absorb_tri_o+=(cur_I-tmp_I);
								tmpobj->I_absorb_o+=(cur_I-tmp_I);
							}
							cur_I=tmp_I;
						}
						if(cur_I<discardI)
						{
							flag_recal=0;
							I_discard=I_discard+cur_I;
						}
					}
					else
					{
						// do nothing since cur_I<discard_I. continue next while-loop and scanf next line.
					}
				}
			}
			fclose(fp);
		}
	}
		

	/*****************************************************************
	*       divide leaf into thin layers, and output ab profile
	*****************************************************************/
	printf("Of all %d files, %d threads actually failed\n",max_file_idx,count_fail);
	sum_p();
	
	int k,num_layers,idx_layer;
	num_layers=atoi(*++argv);
	double sum_I_absorb_chl, tmp_z, opt_ab_profile[100];
	FILE *fout_tri, *fout_srf;
	double sum_I_absorb_all=0.0;

	if((fout_tri=fopen(*++argv,"wb"))==NULL)
	{
		printf("error on open file result tri\n");
	}
	if((fout_srf=fopen(*++argv,"wb"))==NULL)
	{
		printf("error on open file result surface\n");
	}

	//initialize opt_ab_profile
	for(i=0;i<num_layers;i++)
	{
		opt_ab_profile[i]=0.0;
	}

	//summarize light absorption by chloroplasts
	for(i=1;i<=msall_num;i++)
	{
		curobj=p_cell_ms[i-1];
		sum_I_absorb_all+=(curobj->I_absorb_o+curobj->I_absorb_i);
		for(k=0;k<curobj->S->tnum;k++)
		{
			fprintf(fout_tri,"%d %d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,(curobj->S->tri+k)->count_hit_tri,
									(curobj->S->tri+k)->I_absorb_tri_o,(curobj->S->tri+k)->I_absorb_tri_i);
		}
		fprintf(fout_srf,"%d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb_o,curobj->I_absorb_i);
		
		for(j=1;j<=ms_chl_num[i-1];j++)
		{
			curobj=p_chl_ms[i-1][j-1];
			sum_I_absorb_all+=(curobj->I_absorb_o+curobj->I_absorb_i);
			for(k=0;k<curobj->S->tnum;k++)
			{
				fprintf(fout_tri,"%d %d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,
										(curobj->S->tri+k)->count_hit_tri,(curobj->S->tri+k)->I_absorb_tri_o,(curobj->S->tri+k)->I_absorb_tri_i);

				sum_I_absorb_chl+=(curobj->S->tri+k)->I_absorb_tri_i;
				
				tmp_z=((curobj->S->vertex+(curobj->S->tri+k)->n[0])->coordinate[3-1]
					+(curobj->S->vertex+(curobj->S->tri+k)->n[1])->coordinate[3-1]
					+(curobj->S->vertex+(curobj->S->tri+k)->n[2])->coordinate[3-1])/3;
				idx_layer=(int)floor((tmp_z-zmin)/((zmax-zmin)/num_layers));//range 0 ~ num_layers-1
			
				if(idx_layer<0 || idx_layer>=num_layers)
				{
					printf("[Error] z is out of range [zmin, zmax].\n");
					printf("    tmp_z=%e\n",tmp_z);
					printf("    zmin=%e, zmax=%e\n",zmin, zmax);
					exit(1);
				}
				else
				{
					opt_ab_profile[idx_layer]=opt_ab_profile[idx_layer]+(curobj->S->tri+k)->I_absorb_tri_i;
				}
			}
			fprintf(fout_srf,"%d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb_o,curobj->I_absorb_i);
		}

		curobj=p_vac_ms[i-1];
		sum_I_absorb_all+=(curobj->I_absorb_o+curobj->I_absorb_i);
		for(k=0;k<curobj->S->tnum;k++)
		{
			fprintf(fout_tri,"%d %d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,
									(curobj->S->tri+k)->count_hit_tri,(curobj->S->tri+k)->I_absorb_tri_o,(curobj->S->tri+k)->I_absorb_tri_i);
		}
		fprintf(fout_srf,"%d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb_o,curobj->I_absorb_i);
	}

	//////////// nonMS ///////////////////
	for(i=1;i<=nonms_num;i++)
	{
		curobj=p_cell_ns[i-1];
		sum_I_absorb_all+=(curobj->I_absorb_o+curobj->I_absorb_i);
		for(k=0;k<curobj->S->tnum;k++)
		{
			fprintf(fout_tri,"%d %d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,k,
									(curobj->S->tri+k)->count_hit_tri,(curobj->S->tri+k)->I_absorb_tri_o,(curobj->S->tri+k)->I_absorb_tri_i);
		}
		fprintf(fout_srf,"%d %d %d %d %d %e %e\n",curobj->celltype,curobj->cellname,curobj->chlname,curobj->vacuolename,curobj->count_hit,curobj->I_absorb_o,curobj->I_absorb_i);	
	}

	fclose(fout_tri);
	fclose(fout_srf);

	//write ab_profile to a txt file
	FILE *fout_abprofile;
	if((fout_abprofile=fopen(*++argv,"wb"))==NULL)
	{
		printf("error on open file - ab profile\n");
	}
	for(i=0;i<num_layers;i++)
	{
		fprintf(fout_abprofile,"%e ",opt_ab_profile[i]);
	}
	fclose(fout_abprofile);

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
	printf("RT recalculation finished.\n");
	printf("time:%f sec\n",difftime(end,start));
	fprintf(fout_sum,"Time %f\n",difftime(end,start));
	fclose(fout_sum);

	return 1;
}
