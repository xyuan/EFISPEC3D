#include "metis.h"
#include "float.h"

// should be defined with -D in a target of makefile
#define VERBOSE

// debug print
#define CHECK printf("<<< OK, func %s, line %d >>>\n", __FUNCTION__, __LINE__); fflush(stdout);


#define debugmode 0
#define DBG if (debugmode) {printf("file %s, function %s, line %d\n", __FILE__, __FUNCTION__, __LINE__);}
#define DBGR if (debugmode) {printf("rank : %d, file %s, function %s, line %d\n", rank, __FILE__, __FUNCTION__, __LINE__);}

// macros
#define MALLOC(TYPE, NUM) (TYPE*)malloc((NUM)*sizeof(TYPE))
#define CALLOC(TYPE, NUM) (TYPE*)calloc((NUM)*sizeof(TYPE))

//#define STOP(msg) {printf("\n!!! %s -> stop\n\n", msg); fflush(stdout); exit(1);}
#define STOP(msg) {printf("\n!!! %s -> stop\n\n", msg); fflush(stdout);}
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
#define HEXA_8_HDR  "TYPE=C3D8R, ELSET=l"
#define HEXA_27_HDR "TYPE=C3D27R, ELSET=l"
#define QUAD_P_HDR  "TYPE=S4R, ELSET=p"
#define QUAD_F_HDR  "TYPE=S4R, ELSET=f"

// TODO : define weights if no snapshots
#ifdef NO_SNAPSHOTS
#define HEXA_WEIGHT     1
#define QUAD_P_WEIGHT   0
#define QUAD_F_WEIGHT   0
#else
//#define HEXA_WEIGHT 23
//#define QUAD_P_WEIGHT 1
//#define QUAD_F_WEIGHT 48
#define HEXA_WEIGHT     1
#define QUAD_P_WEIGHT   0
#define QUAD_F_WEIGHT   0
#endif

// To use efispec3D numerotation
#define NFACE    6
#define NEDGE   12
#define NCORNER  8

#define FACE1 0
#define FACE2 1
#define FACE3 2
#define FACE4 3
#define FACE5 4
#define FACE6 5

#define EDGE1   0
#define EDGE2   1
#define EDGE3   2
#define EDGE4   3
#define EDGE5   4
#define EDGE6   5
#define EDGE7   6
#define EDGE8   7
#define EDGE9   8
#define EDGE10  9
#define EDGE11 10
#define EDGE12 11

// efispec3D to cubit
#define CORNER1 3
#define CORNER2 2
#define CORNER3 1
#define CORNER4 0
#define CORNER5 7
#define CORNER6 6
#define CORNER7 5
#define CORNER8 4

int edge2efispec[12]       = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
int face2efispec[6]        = {1, 2, 3, 4, 5, 6};
int corner2efispec[8]      = {4, 3, 2, 1, 8, 7, 6, 5};
int cornerefispec2cubit[8] = {3, 2, 1, 0, 7, 6, 5, 4};

// Initial state : {Horiz vect, Vert vect} using nodes 1 or 7 in bottom left corner for each face.
// H or V give direction, + or - give sens : toward up & right is +, - for down and left
// 1 is for i, 2 for j, 3 for k
// ie for face 1, considering it with node1 in the bottom left corner, K is the horizontal vetor directed toward right, j is vertical directed toward up
#define OR_F1 { 3,  2}
#define OR_F2 { 1,  3}
#define OR_F3 {-1, -2}
#define OR_F4 {-3, -1}
#define OR_F5 { 2,  1}
#define OR_F6 {-2, -3}
char face_orientation[6][2] = {OR_F1, OR_F2, OR_F3, OR_F4, OR_F5, OR_F6};

#define F1F1 {1,0,3,2}
#define F1F2 {6,2,3,7}
#define F1F3 {2,6,5,1}
#define F1F4 {0,1,5,4}
#define F1F5 {4,7,3,0}
#define F1F6 {7,4,5,6}

#define F2F1 {1,0,3,2}
#define F2F2 {6,2,3,7}
#define F2F3 {2,6,5,1}
#define F2F4 {0,1,5,4}
#define F2F5 {4,7,3,0}
#define F2F6 {7,4,5,6}

#define F3F1 {3,2,1,0}
#define F3F2 {3,7,6,2}
#define F3F3 {5,1,2,6}
#define F3F4 {5,4,0,1}
#define F3F5 {3,0,4,7}
#define F3F6 {5,6,7,4}

#define F4F1 {0,3,2,1}
#define F4F2 {2,3,7,6}
#define F4F3 {6,5,1,2}
#define F4F4 {1,5,4,0}
#define F4F5 {7,3,0,4}
#define F4F6 {4,5,6,7}

#define F5F1 {3,2,1,0}
#define F5F2 {3,7,6,2}
#define F5F3 {5,1,2,6}
#define F5F4 {5,4,0,1}
#define F5F5 {3,0,4,7}
#define F5F6 {5,6,7,4}

#define F6F1 {1,0,3,2}
#define F6F2 {6,2,3,7}
#define F6F3 {2,6,5,1}
#define F6F4 {0,1,5,4}
#define F6F5 {4,7,3,0}
#define F6F6 {7,4,5,6}

// selon la face, indice du premier noeud trouvé dans l'ordre de parcours cubit -> on regarde le noeud d'en face pour en déduire l'orientation
// ne sert pas dans le code, juste un pense bete pour verifier le code spaghetti ci-dessus
// char node2test[6] = {0, 2, 1, 0, 0, 4}

// dim 1 : num face source
// dim 2 : num face target
// dim 3 : num orientation
// content : first target node corresponding to first source node
// ie:
// for (i=0; i<4; i++)
//    if (connex_nodes[4] == face2faceRot[num_conn_source][num_conn_target][i]) return i;
int face2faceRot[6][6][4]= {  {F1F1, F1F2, F1F3, F1F4, F1F5, F1F6 }, // face 1 source
                              {F2F1, F2F2, F2F3, F2F4, F2F5, F2F6 }, // face 2 source
                              {F3F1, F3F2, F3F3, F3F4, F3F5, F3F6 }, // face 3 source
                              {F4F1, F4F2, F4F3, F4F4, F4F5, F4F6 }, // face 4 source
                              {F5F1, F5F2, F5F3, F5F4, F5F5, F5F6 }, // face 5 source
                              {F6F1, F6F2, F6F3, F6F4, F6F5, F6F6 }  // face 6 source
                           };

// define arbitrary order for points composing an edge to define a sens
// l'arete est dirigée dans le sens du vecteur unité (du repere orthonormé) qui lui est colinéaire :i, j ou k
// attention, les numeros suivants sont les numeros cubits commençants en 0
#define EDGE1_ORDER  {3,2}
#define EDGE2_ORDER  {2,1}
#define EDGE3_ORDER  {0,1}
#define EDGE4_ORDER  {3,0}
#define EDGE5_ORDER  {7,6}
#define EDGE6_ORDER  {6,5}
#define EDGE7_ORDER  {4,5}
#define EDGE8_ORDER  {7,4}
#define EDGE9_ORDER  {3,7}
#define EDGE10_ORDER {2,6}
#define EDGE11_ORDER {1,5}
#define EDGE12_ORDER {0,4}

int edgeOrder[12][2] = {   EDGE1_ORDER,
                           EDGE2_ORDER,
                           EDGE3_ORDER,
                           EDGE4_ORDER,
                           EDGE5_ORDER,
                           EDGE6_ORDER,
                           EDGE7_ORDER,
                           EDGE8_ORDER,
                           EDGE9_ORDER,
                           EDGE10_ORDER,
                           EDGE11_ORDER,
                           EDGE12_ORDER
                       };

typedef enum  {
	NONE   = 0,
	HEXA   = 1,
	QUAD_P = 2,
	QUAD_F = 3,
} elt_t;

typedef enum  {
	FACE   = 0,
	EDGE   = 1,
	CORNER = 2,
} topo_t;

typedef struct {
  real_t x;   
  real_t y;
  real_t z;
} vec3_t  ;

// mesh struct
typedef struct {
  int num_nodes_hexa;          // number of geometric nodes per hexahedron (8 or 27)
  int num_node_per_dim;        // number of geometric nodes per direction (2 or 3)
  int npart;                   // number of procs (parts)
  idx_t ne, nn;                // The # of elements and nodes in the mesh
  idx_t nh, nq_parax, nq_surf; // number of hexa, quad parax, quad fsurf
  idx_t ncon;                  // The number of element balancing constraints (element weights)

  idx_t  *eptr, *eind;         // The CSR-structure storing the nodes in the elements elt -> node
  idx_t  *xadj, *adjncy;       // Elements adjacency graph (CSR format)
  idx_t *xadj_hex, *adjncy_hex;// Hexa only adjacency graph (CSR format)
  vec3_t *ncoords;             // The xyz coordinates of geometric nodes of the elements
  idx_t  *part;                // Hexa partition table
  int	 *layer;               // Layer of elmnts	 
  idx_t  *vwgt;                // The weights of the vertices (hexa elts) of hexa adj graph
  elt_t *types;                // elements types
  float* xcoord;
  float* ycoord;
  float* zcoord;
} mesh_t;

// connection with neighbor information
typedef struct ci_t{
        idx_t num_neigh;// neigbor global number
        int num_conn;   // numface, numedge, numcorner
        int orientation;// connection orientation  0  0  s  v  v  s  v  v
                        // hor  vert	
                        // vv: 01=i, 10=j, 11=k s: 1=+, 0=-
        int defined;    // 1 if this struct is filled
        elt_t type;     // HEXA, QUAD_P, QUAD_F
        struct ci_t* next;
} conn_info;

// element information
typedef struct {
        int globalnum;        // global number
        int localnum;         // local number
        char outer;           // element is outer
	conn_info faces  [ 6];// face neighborhood
	conn_info edges  [12];// edge neighborhood
	conn_info corners[ 8];// corner neighborhood
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

// function defined outside in FORTRAN code : 
void init_gll_number_(int* ihexa,int* ngll_total);

void propagate_gll_nodes_face_(int* ihexa, int* iface, int* ielt_num, int* ielt_face, int* ielt_coty);

void propagate_gll_nodes_edge_(int* ihexa, int* iedge, int* ielt_num, int* ielt_edge, int* ielt_coty);

void propagate_gll_nodes_corner_(int* ihexa, int* inode, int* ielt_num, int* ielt_corner);

void propagate_gll_nodes_quad_(int* ihexa, int* iface, int* iquad, int* global_gll_of_quad, int* nb_quad_p);
