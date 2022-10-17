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

int unique_dis(double *all_dis, int vec_size)
{
	// input: a double array
	// output: number of unique dis values
	int i,j,count_unique=0;
	int *flag;
	flag=calloc(vec_size, sizeof(int));//record non-unique element
	for(i=0;i<vec_size;i++)
	{
		*(flag+i)=0;//initialize
	}
	if(vec_size==1)
	{
		count_unique=1;
	}
	else
	{
		for(i=0;i<vec_size-1;i++)
		{
			if(*(flag+i)==0)//skip detected non-unique element
			{
				count_unique=count_unique+1;
				for(j=i+1;j<vec_size;j++)
				{
					if(*(flag+j)==0)//skip detected non-unique element
					{
						if(fabs(*(all_dis+i)-*(all_dis+j))<DBL_EPSILON)
						{
							*(flag+i)=1;
							*(flag+j)=1;
						}
					}
				}
			}
		}
	}
	free(flag);
	return count_unique;
}

void sum_p()
{
	double sum_I_absorb_all=0.0;
	int i,j;
	Object *curobj;

	//ms
	for(i=1;i<=msall_num;i++)
	{
		curobj=p_cell_ms[i-1];
		sum_I_absorb_all+=curobj->I_absorb_o;
		sum_I_absorb_all+=curobj->I_absorb_i;
			
		for(j=1;j<=ms_chl_num[i-1];j++)
		{
			curobj=p_chl_ms[i-1][j-1];
			sum_I_absorb_all+=curobj->I_absorb_o;
			sum_I_absorb_all+=curobj->I_absorb_i;
		}
		
		curobj=p_vac_ms[i-1];
		sum_I_absorb_all+=curobj->I_absorb_o;
		sum_I_absorb_all+=curobj->I_absorb_i;
	}

	//non_ms
	for(i=1;i<=nonms_num;i++)
	{
		curobj=p_cell_ns[i-1];
		sum_I_absorb_all+=curobj->I_absorb_o;
		sum_I_absorb_all+=curobj->I_absorb_i;
	}

	printf("    discardI:%e\n",I_discard);
	printf("    discardI_Reflect:%e\n",I_discard_Rf);
	printf("    discardI_Transmit:%e\n",I_discard_Tr);
	printf("    total absorb I:%e\n",sum_I_absorb_all);
	printf("    sum_I:%e\n\n",sum_I_absorb_all+I_discard+I_discard_Rf+I_discard_Tr);
}

void VecSub(Vec v1,Vec v2,Vec v3)
{
	//Vec v3;
	v3[0]=v2[0]-v1[0];
	v3[1]=v2[1]-v1[1];
	v3[2]=v2[2]-v1[2];
	//return v3;
}

void VecCross(Vec v1,Vec v2,Vec v3)// v1 corss v2
{
	//Vec v3;
	v3[0]=v1[1]*v2[2]-v1[2]*v2[1];
	v3[1]=v1[2]*v2[0]-v1[0]*v2[2];
	v3[2]=v1[0]*v2[1]-v1[1]*v2[0];
	//return v3;
}

void Normalize(Vec v)
{
	double dist;
	dist=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if(!(dist==0.0))
	{
		v[0]=v[0]/dist;
		v[1]=v[1]/dist;
		v[2]=v[2]/dist;
	}
	else
	{
//		printf("Zero-Length Vector cannot be Normalize\n");
		v[0]=0;
		v[1]=0;
		v[2]=0;
		return;
	}
}

void VecCopy(Vec v1,Vec v2)	// copy from v1 to v2
{
	v2[0]=v1[0];
	v2[1]=v1[1];
	v2[2]=v1[2];
}

double VecDot(Vec v1,Vec v2)
{
	double dotvalue;
	dotvalue=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
	return dotvalue;
}

void IsectCopy(Isect *i1,Isect *i2)
{
	i2->inter_tri=i1->inter_tri;
	i2->inter_nor[0]=i1->inter_nor[0];
	i2->inter_nor[1]=i1->inter_nor[1];
	i2->inter_nor[2]=i1->inter_nor[2];
	i2->inter_pt[0]=i1->inter_pt[0];
	i2->inter_pt[1]=i1->inter_pt[1];
	i2->inter_pt[2]=i1->inter_pt[2];
	i2->inter_obj=i1->inter_obj;
	i2->t=i1->t;
}
