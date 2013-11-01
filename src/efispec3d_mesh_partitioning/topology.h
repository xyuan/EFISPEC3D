// To use efispec3D numerotation
#define FACE1 0
#define FACE2 1
#define FACE3 2
#define FACE4 3
#define FACE5 4
#define FACE6 5

#define EDGE1 0
#define EDGE2 1
#define EDGE3 2
#define EDGE4 3
#define EDGE5 4
#define EDGE6 5
#define EDGE7 6
#define EDGE8 7
#define EDGE9 8
#define EDGE10 9
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

int edge2efispec[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
int face2efispec[6] = {1, 2, 3, 4, 5, 6};
int corner2efispec[8] = {4, 3, 2, 1, 8, 7, 6, 5};
int cornerefispec2cubit[8] = {3, 2, 1, 0, 7, 6, 5, 4};

// Initial state : {Horiz vect, Vert vect} using nodes 1 or 7 in bottom left corner for each face.
// H or V give direction, + or - give sens : toward up & right is +, - for down and left
// 1 is for i, 2 for j, 3 for k
// ie for face 1, considering it with node1 in the bottom left corner, K is the horizontal vetor directed toward right, j is vertical directed toward up
#define OR_F1 {3, 2}
#define OR_F2 {1, 3}
#define OR_F3 {-1, -2}
#define OR_F4 {-3, -1}
#define OR_F5 {2, 1}
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
//	for (i=0; i<4; i++)
//		if (connex_nodes[4] == face2faceRot[num_conn_source][num_conn_target][i]) return i;
int face2faceRot[6][6][4]=	{	{F1F1, F1F2, F1F3, F1F4, F1F5, F1F6	},	// face 1 source
								{F2F1, F2F2, F2F3, F2F4, F2F5, F2F6	},	// face 2 source
								{F3F1, F3F2, F3F3, F3F4, F3F5, F3F6	},	// face 3 source
								{F4F1, F4F2, F4F3, F4F4, F4F5, F4F6	},	// face 4 source
								{F5F1, F5F2, F5F3, F5F4, F5F5, F5F6	},	// face 5 source
								{F6F1, F6F2, F6F3, F6F4, F6F5, F6F6	}	// face 6 source
							};

// define arbitrary order for points composing an edge to define a sens
// l'arete est dirigée dans le sens du vecteur unité (du repere orthonormé) qui lui est colinéaire :i, j ou k
// attention, les numeros suivants sont les numeros cubits commençants en 0
#define EDGE1_ORDER {3,2}
#define EDGE2_ORDER {2,1}
#define EDGE3_ORDER {0,1}
#define EDGE4_ORDER {3,0}
#define EDGE5_ORDER {7,6}
#define EDGE6_ORDER {6,5}
#define EDGE7_ORDER {4,5}
#define EDGE8_ORDER {7,4}
#define EDGE9_ORDER {3,7}
#define EDGE10_ORDER {2,6}
#define EDGE11_ORDER {1,5}
#define EDGE12_ORDER {0,4}

int edgeOrder[12][2] = {	EDGE1_ORDER,
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
							EDGE12_ORDER};



	


