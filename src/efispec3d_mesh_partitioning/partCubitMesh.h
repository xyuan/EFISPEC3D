#include "metis.h"

// should be defined with -D in a target of makefile
#define VERBOSE

// macros
#define MALLOC(TYPE, NUM) (TYPE*)malloc((NUM)*sizeof(TYPE))
#define STOP(msg) fprintf(stderr, "\n!!! %s -> stop\n\n",msg); exit(1);
#ifdef VERBOSE
#define MSG(msg) printf("%s\n",msg);
#else
#define MSG(msg) ;
#endif

#define N_BLOCK_MAX 99
#define LENMAX 500
#define MAX_FILENAME_LEN 70

// might be wrong for 27 nodes elts
#define N_NODES_QUAD 4

// patterns for element type recognition in cubit mesh file
#define HEXA_8_HDR "TYPE=C3D8R, ELSET=l"
#define HEXA_27_HDR "TYPE=C3D27R, ELSET=l"
#define QUAD_P_HDR "TYPE=S4R, ELSET=p"
#define QUAD_F_HDR "TYPE=S4R, ELSET=f"

// TODO : define weights if no snapshots
#ifdef NO_SNAPSHOTS
#define HEXA_WEIGHT 1
#define QUAD_P_WEIGHT 0
#define QUAD_F_WEIGHT 0
#else
//#define HEXA_WEIGHT 23
//#define QUAD_P_WEIGHT 1
//#define QUAD_F_WEIGHT 48
#define HEXA_WEIGHT 1
#define QUAD_P_WEIGHT 0
#define QUAD_F_WEIGHT 0
#endif

typedef enum  {
	NONE = 0,
	HEXA = 1,
	QUAD_P = 2,
	QUAD_F = 3,
} elt_t;

typedef enum  {
	FACE = 0,
	EDGE = 1,
	CORNER = 2,
} topo_t;

typedef struct {
  real_t x;   
  real_t y;
  real_t z;
} vec3_t;

// mesh struct
typedef struct {
  int num_nodes_hexa;			// number of geometric nodes per hexahedron (8 or 27)
  int num_node_per_dim;			// number of geometric nodes per direction (2 or 3)
  int npart;					// number of procs (parts)
  idx_t ne, nn;			        // The # of elements and nodes in the mesh
  idx_t nh, nq_parax, nq_surf;	// number of hexa, quad parax, quad fsurf
  idx_t ncon;           		// The number of element balancing constraints (element weights)

  idx_t  *eptr, *eind;   		// The CSR-structure storing the nodes in the elements elt -> node
  idx_t  *xadj, *adjncy; 		// Elements adjacency graph (CSR format)
  idx_t *xadj_hex, *adjncy_hex; // Hexa only adjacency graph (CSR format)
  vec3_t *ncoords;  			// The xyz coordinates of geometric nodes of the elements
  idx_t  *part;		 			// Hexa partition table
  int	 *layer;		// Layer of elmnts	 
  idx_t  *vwgt;          		// The weights of the vertices (hexa elts) of hexa adj graph
  elt_t *types;					// elements types
  char* prefix;
} mesh_t;

// connection with neighbor information
typedef struct ci_t{
	idx_t num_neigh;// neigbor global number
	int num_conn;	// numface, numedge, numcorner
	int orientation;	// connection orientation			0	0	s	v	v	s	v 	v
						//											hor			vert	
						//									vv: 01=i, 10=j, 11=k	s: 1=+, 0=-
	int defined;		// 1 if this struct is filled
	elt_t type;		// HEXA, QUAD_P, QUAD_F
	struct ci_t* next;
} conn_info;

// element information
typedef struct {
	int globalnum;			// global number
	int localnum;			// local number
	char outer;				// element is outer
	conn_info faces[6]; 	// face neighborhood
	conn_info edges[12];	// edge neighborhood
	conn_info corners[8];	// corner neighborhood
} elt_info_t;

// proc information
typedef struct {
	elt_info_t** local_elts;	
	int nb_elt;
	int nb_ext;
	int nb_quad_p;
	int nb_quad_f;
	int nb_conn;
	char* connex;
	char filename[MAX_FILENAME_LEN];
} proc_info_t;

// convenience struct which hold all information
typedef struct {
	proc_info_t* proc;
	elt_info_t* elt;
} info_t;

