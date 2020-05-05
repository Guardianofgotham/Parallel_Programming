#include <igraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "MaxPerm.h"

void create_graph(igraph_t * g)
{
	int i, nVertices=0, nEdges=0, a, b, wt;
	
	scanf("%d %d",&nVertices, &nEdges);
	printf("Nodes: %d Edges: %d\n",nVertices,nEdges);
	igraph_vector_t v;
	igraph_vector_init (&v, 2*nEdges);
	for(i=0;i<2*nEdges;i+=2){
		scanf("%d %d",&a,&b);
		if(a>=nVertices || b>=nVertices){
			printf("vertex id overflow : %d-%d\n",a,b);
			break;
		}
		VECTOR(v)[i]=a;
		VECTOR(v)[i+1]=b;
	}

	igraph_create(g, &v, 0, 0);
	igraph_simplify(g, 1, 0, 0);
}

void getMembershipFromFile(char filename[], igraph_vector_t * membership)
{
	char line[10000];
	FILE * fp;
	fp = fopen(filename, "r");
	VECTOR(*membership)[0] = 0;
	while(fgets(line, 10000, fp))
	{
		char * pch = strtok(line, " \t\r\n");
		int vertex = atoi(pch);
		pch = strtok(NULL, " \t\r\n");
		int comm = atoi(pch);
		VECTOR(*membership)[vertex] = comm;
	}
	fclose(fp);
}

int main(int argc, char * argv[])
{
	int i, nVertices, maxIt;
	double Netw_perm, result;
	char line[10000];
	FILE * fp, * out;
	igraph_t g;

	out = fopen("output.txt", "w");

	igraph_vector_t membership1, membership2, membership3;
	
	create_graph(&g);
	nVertices = igraph_vcount(&g);
	
	igraph_vector_init(&membership1, nVertices);
	
	maxIt = 10;//Max Iteration
	time_t t = time(NULL);
	Netw_perm=igraph_community_MaxPerm(&g, &membership1, maxIt);
	t = time(NULL) - t;
	printf("Our MaxPerm Implementation took %ld seconds\n",t);
	fprintf(out, "Our MaxPerm Implementation took %ld seconds\n",t);
	fprintf(out, "Network Permanence is %lf\n",Netw_perm );	
	for(i=0;i<nVertices;++i)
	{
		fprintf(out, "%d\t%ld\n",i,(long)VECTOR(membership1)[i]);
	}
	
	
	/*
	To find the NMI with ground truth file:community.dat

	//Compare with ground truth of LFR
	igraph_vector_init(&membership2, nVertices);
	getMembershipFromFile("community.dat", &membership2);	
	igraph_compare_communities(&membership1, &membership2, &result, IGRAPH_COMMCMP_NMI);
	printf("NMI of our MaxPerm Implementation with Ground Truth of LFR is %lf\n", result);

	*/
	
	//Free allocated resources
	igraph_vector_destroy(&membership1);
	igraph_vector_destroy(&membership2);
	
	igraph_destroy(&g);

	return 0;
}
