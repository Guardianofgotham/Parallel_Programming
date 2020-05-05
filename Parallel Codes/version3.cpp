#include "P_MaxPerm.h"
#include "hclib_cpp.h"
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include <inttypes.h>

using namespace hclib;
//getconf LEVEL1_DCACHE_LINESIZE
long long int CountPermananceCalled = 0;
double TotalDurationSpentOnPerm = 0;
long long int permInCurrentIter = 0;

double PermComm(igraph_t * g, igraph_vector_t *neig_v, int * vertex_community, int myComm)
{
	CountPermananceCalled++;
	permInCurrentIter++;
	time_t start = time(NULL);
	int i,j,neighbor,neighbor_2;
	int *E_v,*edge_temp;
	int numerator=0,denominator=0;
	int Etemp_v;

	E_v=(int *) calloc(igraph_vcount(g), sizeof(int));

	edge_temp=(int*) malloc(sizeof(int));

	int I_v=0,Emax_v=0,D_v=0;
	double cin_v=0.0;
	for (i = 0; i < igraph_vector_size (neig_v); ++i)
	{
		neighbor=VECTOR(*neig_v)[i];
		if(vertex_community[neighbor]==myComm)
		{I_v++;}
		else
		{	E_v[vertex_community[neighbor]]++;
			Etemp_v = E_v[vertex_community[neighbor]];
			if(Etemp_v>Emax_v)
				Emax_v=Etemp_v;}
	}
	
	for (i = 0; i < igraph_vector_size (neig_v); ++i)
	{
		neighbor=VECTOR(*neig_v)[i];
		if(vertex_community[neighbor]!=myComm)
			continue;
		for (j = i+1; j < igraph_vector_size (neig_v); ++j)
		{
			neighbor_2=VECTOR(*neig_v)[j];
			if(vertex_community[neighbor_2]!=myComm)
				continue;
			denominator++;
			igraph_get_eid(g, edge_temp, neighbor, neighbor_2, 0, 0);
			if((*edge_temp)!=-1)
			{
				numerator++;
			}
		}

	}
	D_v= i>0?i:1;
	if(Emax_v==0)
		Emax_v=1;

	if(denominator==0)
		cin_v=0.0;
	else
		cin_v=(1.0*numerator)/denominator;

	free(E_v);
	free (edge_temp);
	time_t end = time(NULL);
	TotalDurationSpentOnPerm += end-start;
	return ((1.0*I_v)/(Emax_v*D_v))-1+cin_v;	
}

double PermConditional(igraph_t * g, igraph_vector_t *neig_v, int * vertex_community, int v,int myComm, int changedVertex, int changedVertexCommunity)
{
	
	CountPermananceCalled++;
	permInCurrentIter++;
	time_t start = time(NULL);
	int i,j,neighbor,neighbor_2;
	int *E_v,*edge_temp;
	int numerator=0,denominator=0;
	int Etemp_v;
	
	E_v=(int *) calloc(igraph_vcount(g), sizeof(int));
	
	edge_temp=(int*) malloc(sizeof(int));
	int *new_comms = (int*) calloc(igraph_vector_size (neig_v),sizeof(int)); 
	for(int i =0;i<igraph_vector_size (neig_v);++i){
		int k=VECTOR(*neig_v)[i];
		if(k==changedVertex){
			new_comms[i] = changedVertexCommunity;
			continue;
		}
		new_comms[i] = vertex_community[k];
	}
	
	int I_v=0,Emax_v=0,D_v=0;
	double cin_v=0.0;
	for (i = 0; i < igraph_vector_size (neig_v); ++i)
	{
		neighbor=VECTOR(*neig_v)[i];
		if(new_comms[i]==myComm)
		{
			I_v++;
		}
		else
		{
			E_v[new_comms[i]]++;
			Etemp_v = E_v[new_comms[i]];
			if(Etemp_v>Emax_v)
				Emax_v=Etemp_v;
		}
		if(new_comms[i]!=myComm)
			continue;
	}
	
	for (i = 0; i < igraph_vector_size (neig_v); ++i)
	{
		neighbor=VECTOR(*neig_v)[i];
		if(new_comms[i]!=myComm)
			continue;
		for (j = i+1; j < igraph_vector_size (neig_v); ++j)
		{
			neighbor_2=VECTOR(*neig_v)[j];
			if(new_comms[j]!=myComm)
				continue;
			denominator++;
			igraph_get_eid(g, edge_temp, neighbor, neighbor_2, 0, 0);
			if((*edge_temp)!=-1)
			{
				numerator++;
			}
		}

	}
	
	D_v= i>0?i:1;
	if(Emax_v==0)
		Emax_v=1;

	if(denominator==0)
		cin_v=0.0;
	else
		cin_v=(1.0*numerator)/denominator;
	
	free(new_comms);
	free(E_v);
	free (edge_temp);
	
	time_t end = time(NULL);
	TotalDurationSpentOnPerm += end-start;

	return ((1.0*I_v)/(Emax_v*D_v))-1+cin_v;	
}

double igraph_community_MaxPerm(igraph_t * g, igraph_vector_t * membership, int maxIt)
{
	int i,j,u,v,vertices,Itern,count,prev;
	double Sum,Old_Sum,cur_p,n_p,n_p_neig,Netw_perm;
	double cur_p_neig;
	int *Comm_v;
	int C,Original_C;
	igraph_vector_t *neig;
	vertices = igraph_vcount(g);


	neig = (igraph_vector_t *) malloc(vertices * sizeof(igraph_vector_t));

	for (i = 0; i < vertices; ++i)
	{

		igraph_vector_init(&neig[i], 1);
		igraph_neighbors(g, &neig[i], i, IGRAPH_ALL);
		
	}	
	
	int *vertex_community;
	vertex_community = (int *) malloc(vertices*sizeof(int));	//stores which community a node belongs
	Comm_v = (int *) malloc(vertices*sizeof(int));	//stores which community a node belongs
	for (i = 0; i < vertices; ++i)
	{
		vertex_community[i]=i;
	}
	
	Old_Sum=-INT_MAX;
	Sum = Old_Sum + 1;
	printf("%lf\n", Old_Sum);
	Itern=0;
	while(Sum>Old_Sum && Itern<maxIt)
	{
		Itern++;
		Old_Sum=Sum;
		Sum=0;
		permInCurrentIter=0;
		for (v = 0; v < vertices; ++v)
		{

			cur_p=PermComm(g, &(neig[v]), vertex_community, vertex_community[v]);
			if(cur_p==1)
			{
				Sum+=cur_p;
				continue;
			}
			cur_p_neig=0;
			int * isPresent = (int *)calloc(vertices, sizeof(int));
			int numComV = 0;
			double* curr_Permanence_of_neigs = (double *)calloc(igraph_vector_size (&neig[v]), sizeof(double));
			loop_domain_1d *loop1 = new loop_domain_1d(igraph_vector_size (&neig[v]));
			finish([&](){
				forasync1D(loop1,[&](int k){
					int currNeig = VECTOR(neig[v])[k];
					curr_Permanence_of_neigs[k] = PermComm(g,&(neig[currNeig]),vertex_community,vertex_community[currNeig]);

				},false,FORASYNC_MODE_RECURSIVE);//End of for all Loop
			});//End Of Finish
            delete(loop1);
			for (j = 0; j < igraph_vector_size (&neig[v]); ++j)
			{
				u=VECTOR(neig[v])[j];
				int uCommunity = vertex_community[u];
				cur_p_neig+=curr_Permanence_of_neigs[j];
				if(!isPresent[uCommunity])
				{
					isPresent[uCommunity] = 1;
					Comm_v[numComV++] = uCommunity;
				}
			}
			free(isPresent);
            free(curr_Permanence_of_neigs);
			int myComm = vertex_community[v];
			double* perMananceinAllComm = (double*) calloc(numComV,sizeof(double));
			double* allNeigPermSumWhenCommChanged = (double*) calloc(numComV,sizeof(double));
			loop_domain_1d *loop2 = new loop_domain_1d(numComV);
            loop_domain_1d *loop3 = new loop_domain_1d(numComV);
            //Beggining of Finish
			finish([&](){
				// Begin 1 
				// Calculating Permanence by Moving v in all other neighBouring Communities
                forasync1D(loop2, [&perMananceinAllComm,&g,&vertex_community,&Comm_v,v,&neig](int currIndex){
                    int currComm = Comm_v[currIndex];
                    perMananceinAllComm[currIndex] = PermComm(g,&(neig[v]),vertex_community,currComm);
                },false,FORASYNC_MODE_RECURSIVE);//End of loop2
                // End 1.
				
				// Begin 2
				// Calculating Permanence of neighbours when v is being moved
                
                forasync1D(loop3,[v,&Comm_v,&g,&vertex_community,&allNeigPermSumWhenCommChanged,&neig](int currIndex){
                    int changedVertexCommunity = Comm_v[currIndex];
                    double perms[igraph_vector_size(&(neig[v]))];
                    loop_domain_1d *loop3_inner = new loop_domain_1d(igraph_vector_size(&(neig[v])));
                    finish([&](){
                        forasync1D(loop3_inner,[currIndex,&g,&neig,&perms,&vertex_community,v,changedVertexCommunity](int innerCurrIndex){
                            int u = VECTOR(neig[v])[innerCurrIndex];
                            perms[innerCurrIndex] = PermConditional(g,&(neig[u]),vertex_community,u,vertex_community[u],v,changedVertexCommunity);
                        },false,FORASYNC_MODE_RECURSIVE);// End of loop3_inner
                    });// End of Finish
                    delete(loop3_inner);
                    double sum=0.0;
                    for(int k = 0;k<igraph_vector_size(&(neig[v]));k++){
                        sum+=perms[k];
                    }
                    allNeigPermSumWhenCommChanged[currIndex]=sum;
                },false,FORASYNC_MODE_RECURSIVE); // End of loop 3
               
				
			});//End of Finish
            delete(loop2);
            delete(loop3);
			int maxCom=-1;
			double maxVertexPerm=cur_p;
			double maxNeigPerm=cur_p_neig;
			for(int i=0;i<numComV;i++){
				if(maxVertexPerm+maxNeigPerm <  perMananceinAllComm[i]+allNeigPermSumWhenCommChanged[i]){
					maxVertexPerm=perMananceinAllComm[i];
					maxNeigPerm=allNeigPermSumWhenCommChanged[i];
					maxCom = Comm_v[i];
				}
			}
			if(maxVertexPerm+maxNeigPerm>cur_p+cur_p_neig){
				vertex_community[v]=maxCom;
				cur_p=maxVertexPerm;
				cur_p_neig=maxNeigPerm;
			}
			Sum+=cur_p;
			free(perMananceinAllComm);
			free(allNeigPermSumWhenCommChanged);
		}
		printf("SumQ %lf :: iter %d\n", Sum, Itern);
	}
	Netw_perm=Sum/vertices;
	
	for (i = 0; i < vertices; ++i)
	{
		igraph_vector_destroy(&neig[i]);
	}

	free(neig);

	for (i = 0; i < vertices; ++i)
	{
		VECTOR(*membership)[i] = vertex_community[i];
	}

	free(vertex_community);
    free(Comm_v);
	return Netw_perm;
}
