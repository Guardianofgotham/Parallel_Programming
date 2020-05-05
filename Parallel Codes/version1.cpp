#include "P_MaxPerm.h"
#include "hclib_cpp.h"
using namespace hclib;

double Perm(igraph_t * g, igraph_vector_t *neig_v, int * vertex_community, int v)
{
	int i,j,neighbor,neighbor_2;
	int *E_v,*edge_temp;
	int numerator=0,denominator=0;
	int Etemp_v;

	
	E_v=(int *) calloc(igraph_vcount(g), sizeof(int));

	edge_temp=(int*) malloc(sizeof(int));

	int I_v=0,Emax_v=0,D_v=0;
	double cin_v=0.0;
    for(i = 0; i < igraph_vector_size(neig_v); ++i){
        neighbor=VECTOR(*neig_v)[i];
        // Can be Created in Single Sequential Loop
        // Counting number of Internal Connections
		// Counting Number of Maximum Connection To single Community
        if(vertex_community[neighbor]==vertex_community[v])
		{
			I_v++;
		}
		else
		{
			E_v[vertex_community[neighbor]]++;
			Etemp_v = E_v[vertex_community[neighbor]];
			if(Etemp_v>Emax_v)
				Emax_v=Etemp_v;
            continue;
		}
        // Till here
    }
	for (i = 0; i < igraph_vector_size (neig_v); ++i)
	{
		neighbor=VECTOR(*neig_v)[i];

        if(vertex_community[neighbor]!=vertex_community[v])
            continue;
		for (j = i+1; j < igraph_vector_size (neig_v); ++j)
		{
			neighbor_2=VECTOR(*neig_v)[j];
			if(vertex_community[neighbor_2]!=vertex_community[v])
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

	return ((1.0*I_v)/(Emax_v*D_v))-1+cin_v;	
}


double igraph_community_MaxPerm(igraph_t* g, igraph_vector_t* membership, int maxIt){
    int i,j,u,v,vertices,Itern,count,prev;
	double Sum,Old_Sum,cur_p,n_p,n_p_neig,Netw_perm;
	double cur_p_neig;
	int *Comm_v;
	int C,Original_C;
	igraph_vector_t *neig;
	vertices = igraph_vcount(g);

    neig = (igraph_vector_t *) malloc(vertices * sizeof(igraph_vector_t));

    int *vertex_community;
	vertex_community = (int *) malloc(vertices*sizeof(int));	//stores which community a node belongs
	Comm_v = (int *) malloc(vertices*sizeof(int));	//stores which community a node belongs
	

    loop_domain_1d *loop = new loop_domain_1d (vertices);
    finish([&](){
        forasync1D (loop,[&](int i){

            igraph_vector_init(&neig[i], 1);
		    igraph_neighbors(g, &neig[i], i, IGRAPH_ALL);
            vertex_community[i]=i;
            
            //End of Instructions for all loop
            },false, FORASYNC_MODE_RECURSIVE //End of Arguments of For all
        ); // End of For all

    }); // End of Finish
	
	Old_Sum=-INT_MAX;
	Sum = Old_Sum + 1;
	printf("%lf\n", Old_Sum);
	Itern=0;
	while(Sum>Old_Sum && Itern<maxIt)
	{
		Itern++;
		Old_Sum=Sum;
		Sum=0;
		for (v = 0; v < vertices; ++v)
		{
			// if(v % 1000 == 0)
			// 	printf("vertex = %d\n", v);
			cur_p=Perm(g, &(neig[v]), vertex_community, v);
			if(cur_p==1)
			{
				Sum+=cur_p;
				continue;
			}
			cur_p_neig=0;
			int * isPresent = (int *)calloc(vertices, sizeof(int));
			int numComV = 0;


            // Calculating Permanance of Neigbours Parallely 
            loop_domain_1d *loop = new loop_domain_1d (0,igraph_vector_size (&neig[v]));
            double* perms_of_neig = (double *)calloc(igraph_vector_size(&neig[v]), sizeof(double));
            finish([&](){
                forasync1D(loop,[&](int k){
                    int p =VECTOR(neig[v])[k];
                    perms_of_neig[k] = Perm(g, &(neig[p]), vertex_community, p);
                },false,FORASYNC_MODE_RECURSIVE
                );// End of For all
            }); // End of finish


			for (j = 0; j < igraph_vector_size (&neig[v]); ++j)
			{
				u=VECTOR(neig[v])[j];
				int uCommunity = vertex_community[u];
				cur_p_neig+=perms_of_neig[j];
				if(!isPresent[uCommunity])
				{
					isPresent[uCommunity] = 1;
					Comm_v[numComV++] = uCommunity;
				}
			}
			free(isPresent);
			free(perms_of_neig);
			j=0;
			while(j < numComV)
			{
				C=Comm_v[j];
				j++;
				Original_C=vertex_community[v];

				vertex_community[v]=C;
				n_p=Perm(g, &(neig[v]), vertex_community, v);

				n_p_neig=0;
                
                ///Paraleli Calculate Permanance of neighbours
                loop_domain_1d *loop2 = new loop_domain_1d (igraph_vector_size (&neig[v]));
                double* store_Perms = (double*) calloc(igraph_vector_size (&neig[v]),sizeof(double));
                finish([&](){
                    forasync1D(loop2,[&](int index){
                        int k = VECTOR(neig[v])[index];
                        store_Perms[index] = Perm(g, &(neig[k]), vertex_community, k);
                    }
                    ,false, FORASYNC_MODE_RECURSIVE //End of Arguments for All
                    );//End of For all 
                });//End of Finish
                // O(n^3) is reduced to O(n)
				for (i = 0; i < igraph_vector_size (&neig[v]); ++i)
				{
					u=VECTOR(neig[v])[i];
					
					n_p_neig+=store_Perms[i];

				}
                free(store_Perms);
				if (cur_p + cur_p_neig < n_p + n_p_neig)// Many different kinds of conditions can be put here
				{
					cur_p=n_p;
					cur_p_neig=n_p_neig;	//This line not in original algo
				}
				else
				{
					vertex_community[v]=Original_C;
				}							

			}

			Sum+=cur_p;
		}	
		printf("SumQ %lf :: iter %d\n", Sum, Itern);
	}
	Netw_perm=Sum/vertices;

    loop_domain_1d *loop3 = new loop_domain_1d (vertices);
    finish([&neig, &vertex_community, &membership, &loop3](){
        forasync1D(loop3,[&](int index){
            igraph_vector_destroy(&neig[index]);
            VECTOR(*membership)[index] = vertex_community[index];
        },false, FORASYNC_MODE_RECURSIVE//End of Arguments
        );//End of For all
    });//End of Finish

    free(neig);
	free(vertex_community);
	return Netw_perm;
}