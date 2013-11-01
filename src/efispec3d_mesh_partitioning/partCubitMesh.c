/****************************************************************/
// Read Cubit mesh, partition it with Metis 
// and write cpu mesh files for efispec3D
//
// Does not work with 27 nodes hexahedrons yet (see TODO comments)
// Assumes hexahedrons are written before quads in cubit mesh file
//
// Florent De Martin & David Michea @BRGM 2012
// f.demartin@brgm.fr
// d.michea@brgm.fr
/****************************************************************/

// Standard C include
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// my include 
#include "partCubitMesh.h"
#include "topology.h"

// debug print
#define CHECK printf("<<< OK, func %s, line %d >>>\n", __FUNCTION__, __LINE__); fflush(stdout);

// funct decl
int count_elt(FILE *unit);
void setHexaWeight(mesh_t* mesh);
void setPrefix(char* fileinp, mesh_t* mesh);
void readCubitMesh(char* fileinp, mesh_t* mesh);
info_t* getConnInfo(mesh_t* mesh);
int getProcConn(int ielt, int iproc, proc_info_t* proc_tab, mesh_t* mesh);
void setEltConn(int ielt, int iproc, elt_info_t* elt_tab, mesh_t* mesh, proc_info_t* proc_tab);
void setHexaConnDetail(int elt_num, int neigh_num, int nb_conn, int* connex_nodes, elt_info_t* elt_tab, mesh_t* mesh);
void setQuadConnDetail(int elt_num, int neigh_num, int nb_conn, int* connex_nodes, elt_info_t* elt_tab, mesh_t* mesh);
int findFaceOrientation(int *connex_nodes, int num_conn_source, int num_conn_target);
int findEdgeOrientation(int *connex_nodes, int num_conn_source, int num_conn_target);
void writeProcFiles(info_t* info, mesh_t* mesh);
char *ltrim(char *s);
char *rtrim(char *s);
char *trim(char *s);
void printMemUsage(mesh_t* mesh);
char* rotateVect(char* ortn, char* ret_or, int cw_quarter_rot);
void add_conn(topo_t typecon, elt_info_t* elt_tab, int hexa_src, int numcon_src, idx_t num_neigh, int numcon_neigh, int orientation, elt_t type);
int getnbneigh(topo_t typecon, elt_info_t* elt_tab, int numcon_src);
// TODO
// void freeInfo(info_t* info);
// void freeMesh(mesh_t* mesh);

int main(int argc, char *argv[])
{
	mesh_t mesh; // -> the mesh struct
	idx_t ncommon;
	idx_t pnumflag;
	idx_t edgecuts;
	char* fileinp;

	// Get args
	if(argc == 3) {
		fileinp = argv[1];
		mesh.npart = atoi(argv[2]);
	} else {
		printf("usage : %s mesh.inp nb_part\n", argv[0]);
		exit(1);
	}

	MSG("Read mesh file")
	readCubitMesh(fileinp, &mesh);
	printMemUsage(&mesh);

	MSG("Build complete adjacency graph with 4 nodes connections")
	ncommon = 4;
	pnumflag = 0;
	METIS_MeshToDual(&mesh.ne, &mesh.nn, mesh.eptr, mesh.eind, &ncommon, &pnumflag, &mesh.xadj, &mesh.adjncy);

	MSG("Build hexahedron adjacency graph with 1 node connection")
	ncommon = 1;
	METIS_MeshToDual(&mesh.nh, &mesh.nn, mesh.eptr, mesh.eind, &ncommon, &pnumflag, &mesh.xadj_hex, &mesh.adjncy_hex);

	MSG("Set weight to hexa")
	setHexaWeight(&mesh);
	
	MSG("Part mesh using Metis")
	if (mesh.npart > 1) {
//		METIS_PartGraphKway(&mesh.nh, &mesh.ncon, mesh.xadj_hex, mesh.adjncy_hex, mesh.vwgt, NULL, NULL, &mesh.npart, NULL, NULL, NULL, &edgecuts, mesh.part);
  		METIS_PartGraphRecursive(&mesh.nh, &mesh.ncon, mesh.xadj_hex, mesh.adjncy_hex, mesh.vwgt, NULL, NULL, &mesh.npart, NULL, NULL, NULL, &edgecuts, mesh.part);
	} else {
		memset(mesh.part, 0, mesh.nh*sizeof(idx_t));
	}
		
	MSG("Get complete connectivity information")
	info_t* info = getConnInfo(&mesh);

	MSG("Write mesh files for solver")
	writeProcFiles(info, &mesh);

	// freeInfo(info);
	// freeMesh(&mesh);
	
	MSG("Done.")
	return 0;
}

// write domain files for each proc of efispec3D
void writeProcFiles(info_t* info, mesh_t* mesh)
{
	int iproc, iproc2, ielt, i, j, k;
	proc_info_t* proc_tab = info->proc;
	elt_info_t* elt_tab  = info->elt;
		
	// cf BINARY PROC FILES LAYOUT at the end of file
	for(iproc=0; iproc<mesh->npart; iproc++) {
		FILE* fileout = fopen(proc_tab[iproc].filename, "a");

		fwrite(&proc_tab[iproc].nb_conn, sizeof(int), 1, fileout);									// nb connected proc
		for(iproc2 = 0; iproc2<mesh->npart; iproc2++) {
			if (proc_tab[iproc].connex[iproc2]) fwrite(&iproc2, sizeof(int), 1, fileout);			// liste of connected proc
		}
		
		fwrite(&proc_tab[iproc].nb_elt, sizeof(int), 1, fileout);									// nb hexa (int)
		fwrite(&proc_tab[iproc].nb_ext, sizeof(int), 1, fileout);									// nb outer hexa (int)
		fwrite(&proc_tab[iproc].nb_quad_p, sizeof(int), 1, fileout);								// nb quad p (int)
		fwrite(&proc_tab[iproc].nb_quad_f, sizeof(int), 1, fileout);								// nb quad f (int)
		fwrite(&mesh->num_nodes_hexa, sizeof(int), 1, fileout);										// number of geometric nodes per hexa (int)

		printf("proc %d :\t%d hexa, %d outer\n", iproc, proc_tab[iproc].nb_elt, proc_tab[iproc].nb_ext);
		for (ielt=0; ielt<proc_tab[iproc].nb_elt; ielt++) {											// for each hexa of proc iproc
			elt_info_t* curhex = proc_tab[iproc].local_elts[ielt];
			int globnum = curhex->globalnum;
			// layer
			fwrite(&mesh->layer[globnum], sizeof(int), 1, fileout);						// its layer (uchar)
			// geometry
			for (i=0; i<mesh->num_nodes_hexa; i++) {												// geom coords of hexa geom nodes
				// use efispec node order
				fwrite(&mesh->eind[mesh->eptr[globnum]+corner2efispec[i]-1], sizeof(int), 1, fileout);// numerotation start from 1
			}
			// faces
			for (i=0; i<6; i++) {	
				// 6 faces connectivity
				fwrite(&curhex->faces[i].type, sizeof(int), 1, fileout);							// type of the element connected to face i
				if (curhex->faces[i].type ==  HEXA) {
					int nonZeroNum = elt_tab[curhex->faces[i].num_neigh].localnum + 1;
					fwrite(&nonZeroNum, sizeof(int), 1, fileout);	// local num of neighbor (starting from 1)
					fwrite(&face2efispec[curhex->faces[i].num_conn], sizeof(int), 1, fileout);		// num face of neighbor
					fwrite(&curhex->faces[i].orientation, sizeof(int), 1, fileout);				// orientation
					// DAVID : num_neigh = vieil ordre ???
					fwrite(&mesh->part[curhex->faces[i].num_neigh], sizeof(int), 1, fileout);		// proc of neighbor
				}
			}
			// edges
			for (i=0; i<12; i++) {																	// 12 edges connectivity
				int nb_neighbor = getnbneigh(EDGE, curhex, i);
				fwrite(&nb_neighbor, sizeof(int), 1, fileout);										// nb of neighbors connected to edge i
				
				conn_info* ptr = &curhex->edges[i];
				for (j=0; j<nb_neighbor; j++) {
					fwrite(&ptr->type, sizeof(int), 1, fileout);							// type of the element connected to edge i
					if (ptr->type ==  HEXA) {
						int nonZeroNum = elt_tab[ptr->num_neigh].localnum + 1;
						fwrite(&nonZeroNum, sizeof(int), 1, fileout);	// local num of neighbor
						fwrite(&edge2efispec[ptr->num_conn], sizeof(int), 1, fileout);		// num edge of neighbor
						fwrite(&ptr->orientation, sizeof(int), 1, fileout);				// orientation
						fwrite(&mesh->part[ptr->num_neigh], sizeof(int), 1, fileout);		// proc of neighbor
					} else if (ptr->type != NONE) {
						STOP("edge neighbor should be an hexa")
					}
					ptr = ptr->next;
				}
			}
			// corner
			for (j=0; j<8; j++) {
				// 8 corners connectivity
				i = cornerefispec2cubit[j];															// written in efispec order
				int nb_neighbor = getnbneigh(CORNER, curhex, i);
				fwrite(&nb_neighbor, sizeof(int), 1, fileout);										// nb of neighbors connected to corner i
				
				conn_info* ptr = &curhex->corners[i];
				for (k=0; k<nb_neighbor; k++) {
					fwrite(&ptr->type, sizeof(int), 1, fileout);							// type of the element connected to corner i
					if (ptr->type ==  HEXA) {
						int nonZeroNum = elt_tab[ptr->num_neigh].localnum + 1; 
						fwrite(&nonZeroNum, sizeof(int), 1, fileout);	// local num of neighbor (start from 1)
						fwrite(&corner2efispec[ptr->num_conn], sizeof(int), 1, fileout);					// num corner of neighbor
						// fwrite(&ptr->num_conn, sizeof(int), 1, fileout);					// num corner of neighbor
						fwrite(&mesh->part[ptr->num_neigh], sizeof(int), 1, fileout);		// proc of neighbor
					} else if (ptr->type != NONE) {
						STOP("corner neighbor should be an hexa")
					}
					ptr = ptr->next;
				}
			}
		}
		fclose(fileout);
	}
}

// get all topological informations on all elements. return a filled info_t struct
// write geom nodes in files and overwrite mesh->eind with local node num
info_t* getConnInfo(mesh_t* mesh)
{
	int ielt, iproc, iproc2, inodeptr;
	char filename[100];
	
	// data structure allocation
	proc_info_t* proc_tab = MALLOC(proc_info_t, mesh->npart); // DAVID2FREE
	elt_info_t* elt_tab = MALLOC(elt_info_t, mesh->nh); // DAVID2FREE
	info_t* info = MALLOC(info_t, 1); // DAVID2FREE
	info->proc = proc_tab;
	info->elt = elt_tab;

	// alloc & settings
	for(iproc = 0; iproc<mesh->npart; iproc++) {
		// printf("proc num %d\n", iproc);
		proc_tab[iproc].connex = MALLOC(char, mesh->npart); // DAVID2FREE
		memset(proc_tab[iproc].connex, 0, mesh->npart*sizeof(char));
		proc_tab[iproc].nb_elt = 0;
		proc_tab[iproc].nb_conn = 0;
		proc_tab[iproc].nb_ext = 0;
		proc_tab[iproc].nb_quad_p = 0;
		proc_tab[iproc].nb_quad_f = 0;
		
		sprintf(proc_tab[iproc].filename, "%s.%06d.elements.cpu.%06d.dat", mesh->prefix, mesh->npart, iproc);
	}
	
	// get all connectivity info for each element
	for(ielt=0; ielt<mesh->nh; ielt++) {
		//printf("                           \relt : %d / %d", ielt, mesh->nh);
		iproc = mesh->part[ielt];
		elt_tab[ielt].globalnum = ielt;
		
		// check if elem is outer & set proc connectivity info
		elt_tab[ielt].outer = getProcConn(ielt, iproc, proc_tab, mesh);

		if (elt_tab[ielt].outer) proc_tab[iproc].nb_ext++;

		proc_tab[iproc].nb_elt++;
		// element's neighbors connectivity
		setEltConn(ielt, iproc, elt_tab, mesh, proc_tab);
	}
	
	// local geom nodes
	idx_t* nodeMask;
	nodeMask = MALLOC(idx_t, mesh->nn);
	for(iproc = 0; iproc<mesh->npart; iproc++) {
		memset(nodeMask, 0, mesh->nn*sizeof(idx_t));
		int nb_local_nodes = 0;
		FILE* fileout = fopen(proc_tab[iproc].filename, "w");
		
		// dummy -> overwrite it later (ugly but save memory)
		fwrite(&nb_local_nodes, sizeof(int) ,1 ,fileout);
		fwrite(&nb_local_nodes, sizeof(int) ,1 ,fileout);
		
		for(ielt=0; ielt<mesh->nh; ielt++) {
			if (mesh->part[ielt] == iproc) {
				for(inodeptr = mesh->eptr[ielt]; inodeptr<mesh->eptr[ielt+1]; inodeptr++) {
					
					int globalnodenum = mesh->eind[inodeptr];
					
					if (nodeMask[globalnodenum] == 0) {
						nodeMask[globalnodenum] = ++nb_local_nodes;
						mesh->eind[inodeptr] = nb_local_nodes;	// overwrite eind since we don't use it anymore
						fwrite(&mesh->ncoords[globalnodenum].x, sizeof(real_t) ,1 ,fileout);
						fwrite(&mesh->ncoords[globalnodenum].y, sizeof(real_t) ,1 ,fileout);
						fwrite(&mesh->ncoords[globalnodenum].z, sizeof(real_t) ,1 ,fileout);
					} else {
						int localnumber = nodeMask[globalnodenum];
						mesh->eind[inodeptr] = localnumber;		// overwrite eind since we don't use it anymore
					}
				}	
			}
		} 
		// overwrite number of local nodes
		rewind(fileout);
		int real_t_size = sizeof(real_t);
		fwrite(&real_t_size, sizeof(int) ,1 ,fileout);
		fwrite(&nb_local_nodes, sizeof(int) ,1 ,fileout);
		fclose(fileout);
	}
	free(nodeMask);
	
	// alloc & settings
	int num_outer[mesh->npart];
	int num_inner[mesh->npart];
	for(iproc = 0; iproc<mesh->npart; iproc++) {
		num_outer[iproc] = 0;
		num_inner[iproc] = proc_tab[iproc].nb_ext;
		proc_tab[iproc].local_elts = MALLOC(elt_info_t*, proc_tab[iproc].nb_elt); // DAVID2FREE
	}
	
	// build per proc data struct
	for(ielt=0; ielt<mesh->nh; ielt++) {
		int elt_proc = mesh->part[ielt];
		if (elt_tab[ielt].outer) {
			proc_tab[elt_proc].local_elts[num_outer[elt_proc]] = &elt_tab[ielt];
			elt_tab[ielt].localnum = num_outer[elt_proc]++;
		} else {
			proc_tab[elt_proc].local_elts[num_inner[elt_proc]] = &elt_tab[ielt];
			elt_tab[ielt].localnum = num_inner[elt_proc]++;
		}
	}
	return info;
}

// find topological details on connection between one hexahedral elts and all its neighbors, fill structures
void setEltConn(int ielt, int iproc, elt_info_t* elt_tab, mesh_t* mesh, proc_info_t* proc_tab)
{
	int ineig, inode, inode_n, i;
	int num_node_per_face = mesh->num_node_per_dim*mesh->num_node_per_dim;
	int* connex_nodes;

	connex_nodes = MALLOC(int, (2*num_node_per_face));
	
	for (i=0; i<6; i++) {
		elt_tab[ielt].faces[i].type = NONE;
		elt_tab[ielt].faces[i].next = NULL;
		elt_tab[ielt].faces[i].defined = 0;
	}
	for (i=0; i<12; i++) {
		elt_tab[ielt].edges[i].type = NONE;
		elt_tab[ielt].edges[i].next = NULL;
		elt_tab[ielt].edges[i].defined = 0;
	}
	for (i=0; i<8; i++) {
		elt_tab[ielt].corners[i].type = NONE;
		elt_tab[ielt].corners[i].next = NULL;
		elt_tab[ielt].corners[i].defined = 0;
	}
	
	// loop on hexa neighbors
	for (ineig = mesh->xadj_hex[ielt]; ineig < mesh->xadj_hex[ielt+1]; ineig++) {
		int neigh = mesh->adjncy_hex[ineig];
		int nb_conn = 0;
		
		// get nodes shared by ielt and its neighbor
		for(inode = mesh->eptr[ielt]; inode < mesh->eptr[ielt+1]; inode++) {
			for(inode_n = mesh->eptr[neigh]; inode_n < mesh->eptr[neigh+1]; inode_n++) {
				if (mesh->eind[inode] == mesh->eind[inode_n]) {
					connex_nodes[nb_conn] = inode - mesh->eptr[ielt]; // contact nodes numbers for current elt
					connex_nodes[nb_conn+num_node_per_face] = inode_n - mesh->eptr[neigh]; // associated node numbers for its neighbor
					nb_conn++;
				}
			}
		}
		// get connection topology
		if (nb_conn) {
			setHexaConnDetail(ielt, neigh, nb_conn, connex_nodes, elt_tab, mesh);
		} else {
			STOP("neighbor without contact");
		}
	}
	
	// loop on quad neighbors
	for (ineig = mesh->xadj[ielt]; ineig < mesh->xadj[ielt+1]; ineig++) {
		int neigh = mesh->adjncy[ineig];
		if (neigh >= mesh->nh) {	// quad	
			if (mesh->types[neigh] == QUAD_P) {
				proc_tab[iproc].nb_quad_p++;
			} else {
				proc_tab[iproc].nb_quad_f++;			
			}
			int nb_conn = 0;
			// get nodes shared by ielt and its attached quad
			for(inode = mesh->eptr[ielt]; inode < mesh->eptr[ielt+1]; inode++) {
				for(inode_n = mesh->eptr[neigh]; inode_n < mesh->eptr[neigh+1]; inode_n++) {
					if (mesh->eind[inode] == mesh->eind[inode_n]) {
						connex_nodes[nb_conn] = inode - mesh->eptr[ielt]; // contact node number for current elt
						connex_nodes[nb_conn+num_node_per_face] = inode_n - mesh->eptr[neigh]; // associated node number for its neighbor
						nb_conn++;
					}
				}
			}
			if (nb_conn) {
                                assert(nb_conn == N_NODES_QUAD);
				setQuadConnDetail(ielt, neigh, nb_conn, connex_nodes, elt_tab, mesh);
			} else {
				STOP("neighbor quad without contact");
			}
			// TODO
		}
	}
	free(connex_nodes);
	return;
}

// find topological details on connection between a hexahedral elts and an attached quad, fill structures
void setQuadConnDetail(int elt_num, int neigh_num, int nb_conn, int* connex_nodes, elt_info_t* elt_tab, mesh_t* mesh)
{
	int icase, face, inode, inode2, num_conn_source;
	
	switch(mesh->num_nodes_hexa) {
		case 8 : // 8 geometric nodes hexahedrons
			face = -1;
			for(inode = 0; inode<4; inode++) {
				if (connex_nodes[inode] == CORNER1) {
					for(inode2 = 0; inode2<4; inode2++) {
						switch(connex_nodes[inode2]) {
							case CORNER3 : face = FACE1; goto face_continue;
							case CORNER6 : face = FACE2; goto face_continue;
							case CORNER8 : face = FACE5; goto face_continue;
						}
					}
				} else if(connex_nodes[inode] == CORNER7) {
					for(inode2 = 0; inode2<4; inode2++) {
						switch(connex_nodes[inode2]) {
							case CORNER4 : face = FACE4; goto face_continue;
							case CORNER5 : face = FACE6; goto face_continue;
							case CORNER2 : face = FACE3; goto face_continue;
						}
					}
				}
			}
			face_continue:
			if (face == -1) {
				STOP("error in face search")
			}
			add_conn(FACE, elt_tab, elt_num, face, neigh_num, 0, 0, mesh->types[neigh_num]);
			break;
		case 27 : // 27 geometric nodes hexahedrons : TODO
			STOP("27 nodes elts not implemented yet")
	}
	return;
}

// find topological details on connection between 2 hexahedral elts, fill structures
void setHexaConnDetail(int elt_num, int neigh_num, int nb_conn, int* connex_nodes, elt_info_t* elt_tab, mesh_t* mesh)
{
	int icase, inode, inode2, num_conn_source, num_conn_target;
	
	switch(mesh->num_nodes_hexa) {
		case 8 : // 8 geometric nodes hexahedrons
			switch(nb_conn) {
				case 1 : // 1 node connection (corner)
					// num_conn_source = corner2efispec[connex_nodes[0]];
					// num_conn_target = corner2efispec[connex_nodes[4]];
					num_conn_source = connex_nodes[0];
					num_conn_target = connex_nodes[4];
					add_conn(CORNER, elt_tab, elt_num, num_conn_source, neigh_num, num_conn_target, 0, HEXA);
					break;
				case 2 : // 2 nodes connection (edge)
					for(icase = 0; icase<2; icase++) {
						int edge = -1;
						int offset = icase*4;
						for(inode = 0; inode<2; inode++) {
							switch(connex_nodes[inode+offset]) {
								case CORNER1 : 
									for(inode2 = 0; inode2<2; inode2++) {
										switch(connex_nodes[inode2+offset]) {
											case CORNER2 : 
												edge = EDGE1; goto edge_continue;
											case CORNER4 : 
												edge = EDGE4; goto edge_continue;
											case CORNER5 :
												edge = EDGE9; goto edge_continue;
										}
									}
								case CORNER6 : 
									for(inode2 = 0; inode2<2; inode2++) {
										switch(connex_nodes[inode2+offset]) {
											case CORNER2 : 
												edge = EDGE10; goto edge_continue;
											case CORNER7 : 
												edge = EDGE6; goto edge_continue;
											case CORNER5 : 
												edge = EDGE5; goto edge_continue;
										}
									}
								case CORNER3 : 
									for(inode2 = 0; inode2<2; inode2++) {
										switch(connex_nodes[inode2+offset]) {
											case CORNER2 :
												edge = EDGE2; goto edge_continue;
											case CORNER4 :
												edge = EDGE3; goto edge_continue;
											case CORNER7 : 
												edge = EDGE11; goto edge_continue;
										}
									}
								case CORNER8 : 
									for(inode2 = 0; inode2<2; inode2++) {
										switch(connex_nodes[inode2+offset]) {
											case CORNER7 :
												edge = EDGE7; goto edge_continue;
											case CORNER4 :
												edge = EDGE12; goto edge_continue;
											case CORNER5 :
												edge = EDGE8; goto edge_continue;
										}
									}
							}
						}
						edge_continue:
						if (edge == -1) {
							STOP("error in edge search")
						}
						if (icase==0) {
							num_conn_source = edge;
						} else {
							num_conn_target = edge;
						}
					}
					add_conn(EDGE, elt_tab, elt_num, num_conn_source, neigh_num, num_conn_target, findEdgeOrientation(connex_nodes, num_conn_source, num_conn_target), HEXA);
					break;
				case 4 : // 4 nodes connection (face)
					for(icase = 0; icase<2; icase++) {
						int face = -1;
						int offset = icase*4;
						for(inode = 0; inode<4; inode++) {
							if (connex_nodes[inode+offset] == CORNER1) {
								for(inode2 = 0; inode2<4; inode2++) {
									switch(connex_nodes[inode2+offset]) {
										case CORNER3 : face = FACE1; goto face_continue;
										case CORNER6 : face = FACE2; goto face_continue;
										case CORNER8 : face = FACE5; goto face_continue;
									}
								}
							} else if(connex_nodes[inode+offset] == CORNER7) {
								for(inode2 = 0; inode2<4; inode2++) {
									switch(connex_nodes[inode2+offset]) {
										case CORNER4 : face = FACE4; goto face_continue;
										case CORNER5 : face = FACE6; goto face_continue;
										case CORNER2 : face = FACE3; goto face_continue;
									}
								}
							}
						}
						face_continue:
						if (face == -1) {
							STOP("error in face search")
						}
						if (icase == 0) {
							num_conn_source= face;
						} else {
							num_conn_target= face;
						}
					}
					add_conn(FACE, elt_tab, elt_num, num_conn_source, neigh_num, num_conn_target, findFaceOrientation(connex_nodes, num_conn_source, num_conn_target), HEXA);
					break;
				default : // problem !
					STOP("Error in nb nodes connected")
			}
			break;
		case 27 : // 27 geometric nodes hexahedrons : TODO
			STOP("27 nodes elts not implemented yet")
	}
	return;
}

// return true if element is on the boundary
// as a side effect, build procs adjacency graph
int getProcConn(int ielt, int iproc, proc_info_t* proc_tab, mesh_t* mesh)
{
	int ineig;
	int isOuter = 0;

	for (ineig = mesh->xadj_hex[ielt]; ineig < mesh->xadj_hex[ielt+1]; ineig++) {
		int neigh = mesh->adjncy_hex[ineig];
		int neig_proc = mesh->part[neigh];
		if (neig_proc != iproc) {
			if (!proc_tab[iproc].connex[neig_proc]) {
				proc_tab[iproc].connex[neig_proc] = 1;
				proc_tab[iproc].nb_conn++;
			}
			isOuter = 1;
		}
	}
	return isOuter;
}

// set elts weight for partionning
void setHexaWeight(mesh_t* mesh)
{
	idx_t ihex, ineigh;

	mesh->ncon = 1;
	mesh->vwgt = MALLOC(idx_t, mesh->nh); // DAVID2FREE
	for(ihex=0; ihex<mesh->nh; ihex++) {
		mesh->vwgt[ihex] = HEXA_WEIGHT;
		for (ineigh = mesh->xadj[ihex]; ineigh<mesh->xadj[ihex+1]; ineigh++) {
			idx_t neigh = mesh->adjncy[ineigh];
			if(mesh->types[neigh] == QUAD_P) mesh->vwgt[ihex] += QUAD_P_WEIGHT;
			if(mesh->types[neigh] == QUAD_F) mesh->vwgt[ihex] += QUAD_F_WEIGHT;
		}
	}
	return;
}

void printMemUsage(mesh_t* mesh) {

	double sum = 0.;
	int i=1;
	char* unit;
	const long un = 1;
	
	sum += mesh->npart * sizeof(proc_info_t);
	sum += mesh->nh * sizeof(elt_info_t);
	sum += sizeof(info_t);
	sum += mesh->npart * mesh->npart * sizeof(char);
	sum += mesh->nh * sizeof(elt_info_t*);
	sum += (2*mesh->nh + mesh->ne+1 + mesh->nh * mesh->num_nodes_hexa + (mesh->nq_parax + mesh->nq_surf) * N_NODES_QUAD) * sizeof(idx_t);
	sum += mesh->nh * sizeof(unsigned char);
	sum += mesh->nn * sizeof(vec3_t);
	sum += mesh->ne * sizeof(elt_t);
	sum += mesh->nn * sizeof(idx_t);

	while(sum > un<<(i*10)) i++;
	switch(i) {
		case 1 : unit="Octets"; break;
		case 2 : unit="Kio"; break;
		case 3 : unit="Mio"; break;
		case 4 : unit="Gio"; break;
		case 5 : unit="Tio"; break;
		default: unit="Tio"; i=5;
	}
	printf("This program will use %d %s of memory\n", (int)(sum/(1<<((i-1)*10))),unit);
}

void setPrefix(char* fileinp, mesh_t* mesh)
{
	int offset;
	char *sym, *ptr = fileinp;;
	while(sym = strstr(ptr, "/")) ptr = ++sym;

	// if name is too long, truncate it
	if (strlen(ptr) > (MAX_FILENAME_LEN-30+4)) {
		offset = MAX_FILENAME_LEN-30;
		MSG("prefix too long: truncated")
	} else {
		offset = strlen(ptr)-strlen(".inp");
	}
	ptr[offset] = '\0';
	mesh->prefix = ptr;
}

// read a cubit mesh and fill mesh struct
void readCubitMesh(char* fileinp, mesh_t* mesh)
{
	char line[LENMAX];
	int ielt, inode, itmp;
	int loc_eind = 0;
	int ielt_glo = 0;
	int iblock = 0;
	FILE *unit_inp;
	int number_of_elemnt_block[N_BLOCK_MAX];
	char *hexa8=NULL, *hexa27=NULL;
	char *quadf=NULL, *quadp=NULL;
	int ihexa=0;
	int readquad = 0;

	unit_inp = fopen(fileinp,"r");
	if(unit_inp == NULL)
	{
	  STOP("can't open input file")
	}

	setPrefix(fileinp, mesh);

	//First pass on file *.inp: count number of nodes, hexa, quad_parax and quad_fsurf
	mesh->nh = mesh->nq_parax = mesh->nq_surf = 0;
	while (!feof(unit_inp))
	{
		fgets(line, LENMAX, unit_inp);

		//count number of geometric nodes
		if (strstr(line,"NSET=ALLNODES")) mesh->nn = count_elt(unit_inp);

		//count number of hexa
		if ((hexa8 = strstr(line,HEXA_8_HDR)) || (hexa27 = strstr(line,HEXA_27_HDR))) {
			// set number of node in hexa elts
			mesh->num_nodes_hexa = hexa8?8:27;
			number_of_elemnt_block[iblock] = count_elt(unit_inp);
			mesh->nh += number_of_elemnt_block[iblock++];
			if (readquad) {
				STOP("unexpected : quads are read before hexa in cubit mesh")
			}
		}

		//count number of quad_parax
		if (strstr(line,QUAD_P_HDR)) {
			number_of_elemnt_block[iblock] = count_elt(unit_inp);
			mesh->nq_parax += number_of_elemnt_block[iblock++];
			readquad = 1;
		}

		//count number of quad_fsurf
		if (strstr(line,QUAD_F_HDR)) {
			number_of_elemnt_block[iblock] = count_elt(unit_inp);
			mesh->nq_surf += number_of_elemnt_block[iblock++];
			readquad = 1;
		}

	}
	mesh->ne = mesh->nh + mesh->nq_parax + mesh->nq_surf;
	mesh->ncon = 0;   
	mesh->eptr = MALLOC(idx_t, (mesh->ne+1)); // DAVID2FREE
	mesh->eind = MALLOC(idx_t, (mesh->nh * mesh->num_nodes_hexa + (mesh->nq_parax + mesh->nq_surf) * N_NODES_QUAD)); // DAVID2FREE
	mesh->layer = MALLOC(int, mesh->nh); // DAVID2FREE
	mesh->ncoords = MALLOC(vec3_t, mesh->nn); // DAVID2FREE
	mesh->types = MALLOC(elt_t, mesh->ne); // DAVID2FREE
	mesh->part = MALLOC(idx_t, mesh->nh); // DAVID2FREE
	mesh->num_node_per_dim = (int) cbrt(mesh->num_nodes_hexa);
   
  
    //Second pass on file *.inp: make eind & eptr
	rewind(unit_inp);
	iblock = ihexa = 0;
	while (!feof(unit_inp))
	{
		fgets(line, LENMAX, unit_inp);

		//fill array gnodes with coordinates of geometric nodes
		if (strstr(line,"NSET=ALLNODES"))
			for (inode=0;inode<mesh->nn;inode++)
				fscanf(unit_inp,"%d," "%"SCREAL "," "%"SCREAL "," "%"SCREAL, &itmp, &mesh->ncoords[inode].x, &mesh->ncoords[inode].y, &mesh->ncoords[inode].z); //SCREAL is defined in metis.h

		//fill eptr and eind for hexa
		if ((hexa8 = strstr(line,HEXA_8_HDR)) || (hexa27 = strstr(line,HEXA_27_HDR)))
		{
			// get num of layer
			char* hdr = hexa8?HEXA_8_HDR:HEXA_27_HDR;
			int layer = atoi(rtrim(strstr(line,hdr) + strlen(hdr)));
			for (ielt=0;ielt<number_of_elemnt_block[iblock];ielt++)
			{
				fscanf(unit_inp,"%d",&itmp);
				mesh->eptr[ielt_glo] = loc_eind;
				mesh->types[ielt_glo++] = HEXA;
				mesh->layer[ihexa++] = layer;
				for (inode=0;inode<mesh->num_nodes_hexa;inode++)
				{
					fscanf(unit_inp,",%d",&mesh->eind[loc_eind]);
					mesh->eind[loc_eind++]--;	// because num starts on 1 in cubit meshes
				}
			}
			iblock++;
		}

	    //fill eptr and eind for quad_parax and quad_fsurf
		if ((quadp = strstr(line,QUAD_P_HDR)) || (quadf = strstr(line,QUAD_F_HDR)))
		{
			for (ielt=0;ielt<number_of_elemnt_block[iblock];ielt++)
			{
				fscanf(unit_inp,"%d",&itmp);
				mesh->eptr[ielt_glo] = loc_eind;
				mesh->types[ielt_glo++] = quadp?QUAD_P:QUAD_F;
				for (inode=0;inode<N_NODES_QUAD;inode++)
				{
					fscanf(unit_inp,",%d",&mesh->eind[loc_eind]);
					mesh->eind[loc_eind++]--;	// because num starts on 1 in cubit meshes
				}
			}
			iblock++;
		}
	}
	mesh->eptr[ielt_glo] = loc_eind;
	
	if (mesh->nh == 0) {
		STOP("not a valid mesh file");
	}
   
#ifdef VERBOSE   
   MSG("Mesh information :")
   printf("\t%d nodes\n",mesh->nn);
   printf("\t%d hexa\n",mesh->nh);
   printf("\t%d quad_parax\n",mesh->nq_parax);
   printf("\t%d quad_fsurf\n",mesh->nq_surf);
#endif

	return;
}

char* rotateVect(char* ortn, char* ret_or, int cw_quarter_rot)
{
	// ortn[0] -> horizontal vect
	// ortn[1] -> vertical vect
	// 1=i, 2=j, 3=k
	// sign :	+ : left->right or bottom->up
	//			- : left<-right or bottom<-up
	// ie ortn[] = {1, -3} ->
	//	+-----> i
	//	|
	//	| k
	//	V
	
	char newH, newV, Hsign, Vsign;

	char H = ortn[0];
	char V = ortn[1];
	// -H(00) -> +V(11) -> +H(10) -> -V(01) clockwise rotation. H or V : right digit, sign : left digit
	//    0         3         2         1
	char Hrot = H>0?2:0;
	char Vrot = V>0?3:1;
	int nbrot = cw_quarter_rot%4;	

	Hrot = (4-nbrot+Hrot)&3;
	Hsign = (Hrot>>1)?1:-1;
	Vrot = (4-nbrot+Vrot)&3;
	Vsign = (Vrot>>1)?1:-1;

	
	if (Vrot&1) ret_or[1] = Vsign*abs(V);	// impair -> Vert
	else ret_or[0] = Vsign*abs(V);			// pair -> Hor

	if (Hrot&1) ret_or[1] = Hsign*abs(H);	// impair -> Vert
	else ret_or[0] = Hsign*abs(H);			// pair -> Hor

	return ret_or;
}

int findFaceOrientation(int *connex_nodes, int num_conn_source, int num_conn_target)
{
	char or[2];
	char or_rot[2];
	int i, rot;
	for (i=0; i<4; i++) {
		if (connex_nodes[4] == face2faceRot[num_conn_source][num_conn_target][i]) {
			rot = i;
			break;
		}
	}
	// printf("node s : %d, node t : %d, face2faceRot[%d][%d][%d]=%d\n",connex_nodes[0], connex_nodes[4], num_conn_source+1, num_conn_target+1, i, face2faceRot[num_conn_source][num_conn_target][i]);
	memcpy(or, face_orientation[num_conn_target], 2*sizeof(char));
	
	// horizontal flip (on voit la face depuis l'arriere)
	or[0] = -or[0];
	// -rot : sens inverse parce que vu a l'envers
	rotateVect(or, or_rot, -rot);
	char hsign = or_rot[0]>0?1:0;
	char vsign = or_rot[1]>0?1:0;
	char h=abs(or_rot[0]);
	char v=abs(or_rot[1]);
	char valret = (v&3)|(vsign<<2)|((h&3)<<3)|(hsign<<5);
	// printf("rot : %d, face source %d, face target %d : h:%s%s, v:%s%s, " BYTETOBINARYPATTERN , rot, num_conn_source+1, num_conn_target+1, hsign?"+":"-",(abs(or_rot[0])==1)?"i":(abs(or_rot[0])==2)?"j":"k", vsign?"+":"-", (abs(or_rot[1])==1)?"i":(abs(or_rot[1])==2)?"j":"k", BYTETOBINARY(valret));
	// printf("\n\n");

	return (int) valret;
}

int findEdgeOrientation(int *connex_nodes, int num_edge_source, int num_edge_target)
{	
	int source_order=0;
	if (edgeOrder[num_edge_source][0] == connex_nodes[0]) source_order = 1;
	
	if (edgeOrder[num_edge_target][0] == connex_nodes[4]) {
		if (source_order) return 1; // same direction
	} else {
		if (!source_order) return 1; // same direction
	}
	
	return 0;	// opposite directions
}


void add_conn(topo_t typecon, elt_info_t* elt_tab, int hexa_src, int numcon_src, idx_t num_neigh, int numcon_neigh, int orientation, elt_t type)
{
	conn_info *ptr, *last;

	switch(typecon) {
		case FACE :	if (elt_tab[hexa_src].faces[numcon_src].defined) { // malloc
						fprintf(stderr, "ERROR : face %d of elt %d is already connected\n", numcon_src, hexa_src);
						exit(1);
					} else {
						ptr = &elt_tab[hexa_src].faces[numcon_src];
					}
					break;
		case EDGE : if (elt_tab[hexa_src].edges[numcon_src].defined) { // malloc
						for (ptr = &elt_tab[hexa_src].edges[numcon_src]; ptr; ptr = ptr->next) last = ptr;
						ptr = last->next = MALLOC(conn_info, 1);
					} else {
						ptr = &elt_tab[hexa_src].edges[numcon_src];
					}
					break;
		case CORNER : if (elt_tab[hexa_src].corners[numcon_src].defined) { // malloc
						for (ptr = &elt_tab[hexa_src].corners[numcon_src]; ptr; ptr = ptr->next) last = ptr;
						ptr = last->next = MALLOC(conn_info, 1);
					} else {
						ptr = &elt_tab[hexa_src].corners[numcon_src];
					}
					break;
	}
	ptr->num_neigh = num_neigh;
	ptr->num_conn = numcon_neigh;
	ptr->orientation = orientation;
	ptr->type = type;
	ptr->next = NULL;
	ptr->defined = 1;
	return;
}

int getnbneigh(topo_t typecon, elt_info_t* elt_tab, int numcon_src)
{
	conn_info *ptr;
	int count=0;
	switch(typecon) {
		case FACE :	printf("WARNING : funct getnbneigh() should not be called for faces\n");
					count++;
					break;
		case EDGE : if (elt_tab->edges[numcon_src].defined) {
						for (ptr = &elt_tab->edges[numcon_src]; ptr; ptr = ptr->next) count++;
					}
					break;
		case CORNER : if (elt_tab->corners[numcon_src].defined) {
						for (ptr = &elt_tab->corners[numcon_src]; ptr; ptr = ptr->next) count++;
					}
					break;
	}
	return count;
}

//This function count the number of lines to compute the number of nodes, etc.
int count_elt(FILE *unit)
{
   const char EOL = '\n';
   char c;
   int ne = 0;

   while(1)
   {
      c = getc(unit);
      if (c == EOL)
         ne++;
      else if (c == '*')
         break;   
   }
   return ne; 
}

char *ltrim(char *s)
{
    while(isspace(*s)) s++;
    return s;
}

char *rtrim(char *s)
{
    char* back = s + strlen(s);
    while(isspace(*--back));
    *(back+1) = '\0';
    return s;
}

char *trim(char *s)
{
    return rtrim(ltrim(s)); 
}

// TODO
/* 
void freeInfo(info_t* info)
{
}
void freeMesh(mesh_t* mesh)
{
}
*/

/*********************************************************************************
BINARY PROC FILES LAYOUT : 
(all num of faces, edges & corners use efispec num convention)

Sizeof(real_t)                          (int)
<nb_nodes>
	<node1.x> <node1.y> <node1.z>   (real_t)
	<node2.x> <node2.y> <node2.z>   (real_t)
	...
<nb_connected_proc>                     (int)
	<procA>                         (int)
	<procB>                         (int)
	... (nb_connected_proc entries)
<nb_hexa>                               (int)
<nb_outer_hexa>                         (int)
<nb_quad_p>                             (int)
<nb_quad_f>                             (int)
<nb_geom_node_per_hexa>                 (int)
	<hexa_layer>                    (int)
		// efispec order & indices starting from 1
		node_num0
		node_num1
		... (nb_geom_node_per_hexa entries)
	face0 :
		<neighbor type> (cf elt_t in .h)    (int)
		if (type == HEXA)
			<local num of neighbor>     (int)
			<num  face of neighbor>     (int)
			<orientation>               (int) -> ... 0 0 s v v s v v + front 0 padding because int 32bits
                                                                 HHHHH VVVVV -> HHHH: horizontal vect, VVVV: vertical vect
                                                                 s (sign) : 1="+", 0="-"
                                                                 vv : on 2 bits : direction : 01="i"(1), 10="j"(2), 11="k"(3)
			<proc of neighbor>          (int)
	face1 :
	... (6 faces)
	
	edge0 : 
		nb neighbors for edge0
			<neighbor type> (cf elt_t in .h)   (int)
			if (type == HEXA)
				<local num of neighbor>    (int)
				<num edge of neighbor>     (int)
				<orientation>              (int)
				<proc of neighbor>         (int)
		... nb neighbors times
	edge1 :
	... (12 edges)
	
	corner0 :
		nb neighbors for corner0
			<neighbor type> (cf elt_t in .h)   (int)
			if (type == HEXA)
				<local num of neighbor>    (int)
				<num face of neighbor>     (int)
				<proc of neighbor>         (int)
		... nb neighbors times
	corner1 : 
	... (8 corners)
	
... nb_hexa entries. nb_outer_hexa first ones are for outer elements
**********************************************************************************/	
