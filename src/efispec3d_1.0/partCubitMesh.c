////!=====================================================================================================================================
////!                                      EFISPEC3D                                             !
////!                            (Elements FInis SPECtraux 3D)                                   !
////!                                                                                            !
////!                               http://efispec.free.fr                                       !
////!                                                                                            !
////!                                                                                            !
////!                            This file is part of EFISPEC3D                                  !
////!              Please refer to http://efispec.free.fr if you use it or part of it            !
////!                                                                                            !
////!                                                                                            !
////!1 ---> French License: CeCILL V2                                                            !
////!                                                                                            !
////!         Copyright BRGM 2009 "contributeurs : Florent  DE MARTIN (BRGM)                     !
////!                                              David    MICHEA    (BRGM)                     !
////!                                              Philippe THIERRY   (Intel)"                   !
////!                                                                                            !
////!         Contact: f.demartin at brgm.fr                                                     !
////!                                                                                            !
////!         Ce logiciel est un programme informatique servant a resoudre l'equation du         !
////!         mouvement en trois dimensions via une methode des elements finis spectraux.        !
////!                                                                                            !
////!         Ce logiciel est regi par la licence CeCILL soumise au droit francais et            !
////!         respectant les principes de diffusion des logiciels libres. Vous pouvez            !
////!         utiliser, modifier et/ou redistribuer ce programme sous les conditions de la       !
////!         licence CeCILL telle que diffusee par le CEA, le CNRS et l'INRIA sur le site       !
////!         "http://www.cecill.info".                                                          !
////!                                                                                            !
////!         En contrepartie de l'accessibilite au code source et des droits de copie, de       !
////!         modification et de redistribution accordes par cette licence, il n'est offert      !
////!         aux utilisateurs qu'une garantie limitee. Pour les memes raisons, seule une        !
////!         responsabilite restreinte pese sur l'auteur du programme, le titulaire des         !
////!         droits patrimoniaux et les concedants successifs.                                  !
////!                                                                                            !
////!         A cet egard l'attention de l'utilisateur est attiree sur les risques associes      !
////!         au chargement, a l'utilisation, a la modification et/ou au developpement et a      !
////!         la reproduction du logiciel par l'utilisateur etant donne sa specificite de        !
////!         logiciel libre, qui peut le rendre complexe a manipuler et qui le reserve donc     !
////!         a des developpeurs et des professionnels avertis possedant des connaissances       !
////!         informatiques approfondies. Les utilisateurs sont donc invites a charger et        !
////!         tester l'adequation du logiciel a leurs besoins dans des conditions permettant     !
////!         d'assurer la securite de leurs systemes et ou de leurs donnees et, plus            !
////!         generalement, a l'utiliser et l'exploiter dans les memes conditions de             !
////!         securite.                                                                          !
////!                                                                                            !
////!         Le fait que vous puissiez acceder a cet en-tete signifie que vous avez pris        !
////!         connaissance de la licence CeCILL et que vous en avez accepte les termes.          !
////!                                                                                            !
////!                                                                                            !
////!2 ---> International license: GNU GPL V3                                                    !
////!                                                                                            !
////!         EFISPEC3D is a computer program that solves the three-dimensional equations of     !
////!         motion using a finite spectral-element method.                                     !
////!                                                                                            !
////!         Copyright (C) 2009 Florent DE MARTIN                                               !
////!                                                                                            !
////!         Contact: f.demartin at brgm.fr                                                     !
////!                                                                                            !
////!         This program is free software: you can redistribute it and/or modify it under      !
////!         the terms of the GNU General Public License as published by the Free Software      !
////!         Foundation, either version 3 of the License, or (at your option) any later         !
////!         version.                                                                           !
////!                                                                                            !
////!         This program is distributed in the hope that it will be useful, but WITHOUT ANY    !
////!         WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    !
////!         PARTICULAR PURPOSE. See the GNU General Public License for more details.           !
////!                                                                                            !
////!         You should have received a copy of the GNU General Public License along with       !
////!         this program. If not, see http://www.gnu.org/licenses/.                            !
////!                                                                                            !
////!                                                                                            !
////!3 ---> Third party libraries                                                                !
////!                                                                                            !
////!         EFISPEC3D uses the following of third party libraries:                             !
////!                                                                                            !
////!           --> METIS 5.1.0 under the Apache License Version 2.0 (compatible with GNU GPL V3)!
////!               see http://glaros.dtc.umn.edu/gkhome/metis/metis/overview                    !
////!                                                                                            !
////!           --> Lib_VTK_IO under GNU GPL V3 License                                          !
////!               see S. Zaghi's website: https://github.com/szaghi/Lib_VTK_IO                 !
////!                                                                                            !
////!           --> INTERP_LINEAR under GNU GPL License                                          !
////!               see J. Burkardt website: http://people.sc.fsu.edu/~jburkardt/                !
////!                                                                                            !
////!                                                                                            !
////!4 ---> Related Articles                                                                     !
////!                                                                                            !
////!         De Martin, F., Matsushima, M., Kawase, H. (BSSA, 2013)                             !
////!            Impact of geometric effects on near-surface Green's functions                   !
////!            doi:10.1785/0120130039                                                          !
////!                                                                                            !
////!         Aochi, H., Ducellier, A., Dupros, F., Delatre, M., Ulrich, T., De Martin, F.,      !
////!         Yoshimi, M., (Pure Appl. Geophys. 2013)                                            !
////!            Finite Difference Simulations of Seismic Wave Propagation for the 2007 Mw 6.6   !
////!            Niigata-ken Chuetsu-Oki Earthquake: Validity of Models and Reliable             !
////!            Input Ground Motion in the Near-Field                                           !
////!            doi:10.1007/s00024-011-0429-5                                                   !
////!                                                                                            !
////!         De Martin, F. (BSSA, 2011)                                                         !
////!            Verification of a Spectral-Element Method Code for the Southern California      !
////!            Earthquake Center LOH.3 Viscoelastic Case                                       !
////!            doi:10.1785/0120100305                                                          !
////!                                                                                            !
////!=====================================================================================================================================

/** \file partCubitMesh.c 
*   \brief This file contains C subroutines to partition a mesh with (<a href="http://glaros.dtc.umn.edu/gkhome/metis/metis/overview" target="_blank">METIS-5.1.0 library</a>) 
*          and to initialize GLL nodes global numbering.
*/

// Standard C include
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>

// my include 
#include "partCubitMesh.h"
#include "mpi.h"


/********************************************************** FUNCTIONS DECLARATION **********************************************************/
int count_elt(FILE *unit);

void setHexaWeight(mesh_t *mesh);

void readCubitMesh_1(char *fileinp, mesh_t *mesh, int *number_of_elemnt_block);

void readCubitMesh_2(char *fileinp, mesh_t *mesh, int *number_of_elemnt_block);

void read_part_cubit_mesh_(char  *fileinp
                         , int   *ig_ncpu
                         , int   *ig_myrank
                         , int   *size_real_t
                         , int   *ig_mesh_nnode
                         , int   *ig_ncpu_neighbor
                         , int   *ig_nhexa
                         , int   *ig_nhexa_outer
                         , int   *ig_nquad_parax
                         , int   *ig_nquad_fsurf
                         , int   *ig_hexa_nnode
                         , int   *ig_quad_nnode
                         , float *rg_mesh_xmax
                         , float *rg_mesh_ymax
                         , float *rg_mesh_zmax
                         , float *rg_mesh_xmin
                         , float *rg_mesh_ymin
                         , float *rg_mesh_zmin);

int fill_mesh_arrays_     (int *ig_ncpu
                         , int *ig_myrank
                         , int *ig_ngll_total
                         , int *ig_nneighbor_all_kind
                         , int *cpu_neighbor
                         , int *ig_cpu_neighbor_info
                         , int *ig_hexa_gnode_glonum
                         , int *ig_quadp_gnode_glonum
                         , int *ig_quadf_gnode_glonum
                         , int *ig_hexa_gll_glonum
                         , int *ig_quadp_gll_glonum
                         , int *ig_quadf_gll_glonum
                         , int *ig_quadp_neighbor_hexa
                         , int *ig_quadp_neighbor_hexaface
                         , int *ig_hexa_material_number);

void fill_efispec_arrays( info_t *info
                        , mesh_t *mesh
                        , int     rank
                        , int    *ig_ngll_total
                        , int    *ig_nneighbor_all_kind
                        , int    *cpu_neighbor
                        , int    *ig_cpu_neighbor_info
                        , int    *ig_hexa_gnode_glonum
                        , int    *ig_quadp_gnode_glonum
                        , int    *ig_quadf_gnode_glonum
                        , int    *ig_hexa_gll_glonum
                        , int    *ig_quadp_gll_glonum
                        , int    *ig_quadf_gll_glonum
                        , int    *ig_quadp_neighbor_hexa
                        , int    *ig_quadp_neighbor_hexaface
                        , int    *ig_hexa_material_number);

info_t *getConnInfo(mesh_t *mesh, int rank, int *ig_mesh_nnode);

int getProcConn(int ielt, int iproc, proc_info_t *proc_tab, mesh_t *mesh);

void setEltConn(int ielt, int iproc, elt_info_t *elt_tab, mesh_t *mesh, proc_info_t *proc_tab);

void setHexaConnDetail(int elt_num, int neigh_num, int nb_conn, int *connex_nodes, elt_info_t *elt_tab, mesh_t *mesh);

void setQuadConnDetail(int elt_num, int neigh_num, int nb_conn, int *connex_nodes, elt_info_t *elt_tab, mesh_t *mesh);

int findFaceOrientation(int *connex_nodes, int num_conn_source, int num_conn_target);

int findEdgeOrientation(int *connex_nodes, int num_conn_source, int num_conn_target);

void writeProcFiles(info_t *info, mesh_t *mesh);

char *ltrim(char*s);
char *rtrim(char*s);
char *trim (char*s);

void printMemUsage(mesh_t *mesh);

char* rotateVect(char *ortn, char *ret_or, int cw_quarter_rot);

void add_conn(topo_t typecon, elt_info_t *elt_tab, int hexa_src, int numcon_src, idx_t num_neigh, int numcon_neigh, int orientation, elt_t type);

int getnbneigh(topo_t typecon, elt_info_t *elt_tab, int numcon_src);

void free_all_1();

void free_all_2();

void get_domain_min_max(mesh_t *mesh, float *rg_mesh_xmax, float *rg_mesh_ymax, float *rg_mesh_zmax, float *rg_mesh_xmin, float *rg_mesh_ymin, float *rg_mesh_zmin);

// DGN : UGLY workaround
void copy_coords_(float *fortran_array, int *nbelem, int *dim);

// Fortran routines
void mod_init_mesh_mp_init_gll_number_(int*,int*);

void mod_init_mesh_mp_propagate_gll_nodes_face_(int*, int*, int*, int*, int*);

void mod_init_mesh_mp_propagate_gll_nodes_edge_(int*, int*, int*, int*, int*);

void mod_init_mesh_mp_propagate_gll_nodes_corner_(int*, int*, int*, int*);

void mod_init_mesh_mp_propagate_gll_nodes_quad_(int*, int*, int*, int*, int*);

void flush()
{ 
   fflush(stdout);
   fflush(stderr);
}

/********************************************************** GLOBAL VARIABLES **********************************************************/
static mesh_t *mesh_ptr;
static info_t *info_ptr;

float *xcoord;
float *ycoord;
float *zcoord;

/******************************************************** FUNCTIONS DEFINITION ********************************************************/
void copy_coords_(float* fortran_array, int* nbelem, int* dim)
{ 
   float* ptr;
   switch(*dim)
   {
      case 1 : ptr = xcoord; break; 
      case 2 : ptr = ycoord; break;
      case 3 : ptr = zcoord; break;
   }
   
   memcpy(fortran_array, ptr, *nbelem * sizeof(float));
   
}

void free_all_1()
{ 
   free(mesh_ptr->xadj);
   free(mesh_ptr->adjncy);
   free(mesh_ptr->xadj_hex);
   free(mesh_ptr->adjncy_hex);
   free(mesh_ptr->ncoords);
}

void free_all_2()
{ 
   int nbelt = info_ptr->proc->nb_elt;

   free(info_ptr->proc->local_elts);
   free(info_ptr->proc->connex);
   free(info_ptr->proc);

   for(int i=0; i< nbelt; i++) {
      elt_info_t elt = info_ptr->elt[i];

      for (int j=0; j<NEDGE; j++) {
         // conn_info* ptr = &elt->edges[j];
         conn_info* ptr = &elt.edges[j];
         conn_info* next;
         int first = 1;
         while(ptr)
         {   next = ptr->next;
            if (!first) free(ptr);
            ptr = next;
            first = 0;
         }
      }
      for (int j=0; j<NCORNER; j++) {
         // conn_info* ptr = &elt->corners[j];
         conn_info* ptr = &elt.corners[j];
         conn_info* next;
         int first = 1;
         while(ptr)
         {   next = ptr->next;
            if (!first) free(ptr);
            ptr = next;
            first = 0;
         }
      }
   }

   free(info_ptr->elt);
   free(info_ptr);

   free(mesh_ptr->eptr);
   free(mesh_ptr->eind);
   free(mesh_ptr->part);
   free(mesh_ptr->layer);
   free(mesh_ptr->vwgt);
   free(mesh_ptr->types);

   // if no fortran copy of these arrays yet --> do not free
   free(xcoord);
   free(ycoord);
   free(zcoord);

   free(mesh_ptr);
}

/******************************************************** FUNCTIONS *******************************************************************/
void read_part_cubit_mesh_( char * fileinp
                         , int  * ig_ncpu
                         , int  * ig_myrank
                         , int  * size_real_t
                         , int  * ig_mesh_nnode
                         , int  * ig_ncpu_neighbor
                         , int  * ig_nhexa
                         , int  * ig_nhexa_outer
                         , int  * ig_nquad_parax
                         , int  * ig_nquad_fsurf
                         , int  * ig_hexa_nnode
                         , int  * ig_quad_nnode
                         , float* rg_mesh_xmax
                         , float* rg_mesh_ymax
                         , float* rg_mesh_zmax
                         , float* rg_mesh_xmin
                         , float* rg_mesh_ymin
                         , float* rg_mesh_zmin)
{ 
   int rank;
   mesh_t* mesh; // -> the mesh structure
   idx_t ncommon;
   idx_t pnumflag;
   idx_t edgecuts;
   int size[5];
   int number_of_elemnt_block[N_BLOCK_MAX];

   mesh        = MALLOC(mesh_t, 1);
   mesh->npart = *ig_ncpu;
   rank        = *ig_myrank;


   if (rank == 0) {

      MSG("")
      MSG("Read mesh file")
      readCubitMesh_1(fileinp, mesh, number_of_elemnt_block);

      size[0] = mesh->nn;
      size[1] = mesh->nh;
      size[2] = mesh->nq_parax;
      size[3] = mesh->nq_surf;
      size[4] = mesh->num_nodes_hexa;

      MPI_Bcast( size, 5, MPI_INT, 0, MPI_COMM_WORLD);

   } else {

      MPI_Bcast( size, 5, MPI_INT, 0, MPI_COMM_WORLD);
      
      mesh->nn             = size[0];
      mesh->nh             = size[1];
      mesh->nq_parax       = size[2];
      mesh->nq_surf        = size[3];
      mesh->num_nodes_hexa = size[4];

   }

   // allocation
   mesh->ne               = mesh->nh + mesh->nq_parax + mesh->nq_surf;
   mesh->ncon             = 0;
   mesh->eptr             = MALLOC(idx_t, (mesh->ne+1)); // DAVID2FREE
   mesh->eind             = MALLOC(idx_t, (mesh->nh * mesh->num_nodes_hexa + (mesh->nq_parax + mesh->nq_surf) * N_NODES_QUAD)); // DAVID2FREE
   mesh->layer            = MALLOC(int, mesh->nh);       // DAVID2FREE
   mesh->ncoords          = MALLOC(vec3_t, mesh->nn);    // DAVID2FREE
   mesh->types            = MALLOC(elt_t, mesh->ne);     // DAVID2FREE
   mesh->part             = MALLOC(idx_t, mesh->nh);     // DAVID2FREE
   mesh->num_node_per_dim = (int) cbrt(mesh->num_nodes_hexa);

   // MPI type for vec3_t
   MPI_Datatype Vec3_mpi_t;
   MPI_Type_contiguous(3, MPI_FLOAT, &Vec3_mpi_t);
   MPI_Type_commit(&Vec3_mpi_t);

   if (rank == 0) {

      readCubitMesh_2(fileinp, mesh, number_of_elemnt_block);
//    printMemUsage(mesh);

   }

   MPI_Bcast( mesh->eptr   , mesh->ne+1, MPI_INT   , 0, MPI_COMM_WORLD);
   MPI_Bcast( mesh->eind   , (mesh->nh * mesh->num_nodes_hexa + (mesh->nq_parax + mesh->nq_surf) * N_NODES_QUAD), MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast( mesh->layer  , mesh->nh  , MPI_INT   , 0, MPI_COMM_WORLD);
   MPI_Bcast( mesh->ncoords, mesh->nn  , Vec3_mpi_t, 0, MPI_COMM_WORLD);
   MPI_Bcast( mesh->types  , mesh->ne  , MPI_INT   , 0, MPI_COMM_WORLD);

   get_domain_min_max(  mesh, rg_mesh_xmax, rg_mesh_ymax, rg_mesh_zmax, rg_mesh_xmin, rg_mesh_ymin, rg_mesh_zmin);

   if (rank == 0) MSG("Build complete adjacency graph with 4 nodes connections")
   ncommon  = 4;
   pnumflag = 0;
   METIS_MeshToDual(&mesh->ne, &mesh->nn, mesh->eptr, mesh->eind, &ncommon, &pnumflag, &mesh->xadj, &mesh->adjncy);

   if (rank == 0) MSG("Build hexahedron adjacency graph with 1 node connection")
   ncommon = 1;
   METIS_MeshToDual(&mesh->nh, &mesh->nn, mesh->eptr, mesh->eind, &ncommon, &pnumflag, &mesh->xadj_hex, &mesh->adjncy_hex);

   if (rank == 0) MSG("Set weight to hexa")
   setHexaWeight(mesh);
      
   if (rank == 0) MSG("Part mesh using Metis")

   if (mesh->npart > 1) {

      METIS_PartGraphRecursive(&mesh->nh, &mesh->ncon, mesh->xadj_hex, mesh->adjncy_hex, mesh->vwgt, NULL, NULL, &mesh->npart, NULL, NULL, NULL, &edgecuts, mesh->part);
      int numelem = 0;

      for (int hh=0; hh<mesh->nh; hh++) {
         if (mesh->part[hh] == rank) numelem++;
      }

   } else {

      memset(mesh->part, 0, mesh->nh*sizeof(idx_t));

   }

   // rg_[xyz]_coord_geom_nodes allocated and filled in getConnInfo()
   if (rank == 0) MSG("Get complete connectivity information")

   // allocation coord_geom_nodes -> peut-etre redecouper ici !!!
   info_t* info = getConnInfo(mesh, rank, ig_mesh_nnode);

   //sizes 4 efispec
   *size_real_t                  = sizeof(real_t);
   *ig_nhexa        = info->proc->nb_elt;
   *ig_nhexa_outer  = info->proc->nb_ext;
   *ig_nquad_parax  = info->proc->nb_quad_p;
   *ig_nquad_fsurf  = info->proc->nb_quad_f;
   *ig_hexa_nnode = mesh->num_nodes_hexa;
   *ig_ncpu_neighbor   = info->proc->nb_conn;
   *ig_quad_nnode = mesh->num_node_per_dim*mesh->num_node_per_dim;

   // DAVID : ptr a sauvegarder en fortran
   mesh_ptr = mesh;
   info_ptr = info;

   // on garde un pointeur pour liberer plus tard (car on copie coté fortran, sinon c'est inutile !)
   mesh->xcoord = xcoord;
   mesh->ycoord = ycoord;
   mesh->zcoord = zcoord;

   // free what's not needed anymore to avoid a memory consumption pick as much as possible
   free_all_1();

   return;
}

void get_domain_min_max(mesh_t* mesh, float* rg_mesh_xmax, float* rg_mesh_ymax, float* rg_mesh_zmax, float* rg_mesh_xmin, float* rg_mesh_ymin, float* rg_mesh_zmin)
{
   float xmin, ymin, zmin;
   float xmax, ymax, zmax;
   xmin = ymin = zmin = +FLT_MAX;
   xmax = ymax = zmax = -FLT_MAX;

   for (int node=0; node<mesh->nn; node++) {
         if (mesh->ncoords[node].x < xmin) xmin = mesh->ncoords[node].x;
         if (mesh->ncoords[node].y < ymin) ymin = mesh->ncoords[node].y;
         if (mesh->ncoords[node].z < zmin) zmin = mesh->ncoords[node].z;

         if (mesh->ncoords[node].x > xmax) xmax = mesh->ncoords[node].x;
         if (mesh->ncoords[node].y > ymax) ymax = mesh->ncoords[node].y;
         if (mesh->ncoords[node].z > zmax) zmax = mesh->ncoords[node].z;
   }

   *rg_mesh_xmax = xmax;
   *rg_mesh_ymax = ymax;
   *rg_mesh_zmax = zmax;

   *rg_mesh_xmin = xmin;
   *rg_mesh_ymin = ymin;
   *rg_mesh_zmin = zmin;
}

int fill_mesh_arrays_(int* nb_part_,
                      int* rank_,
                      int* ig_ngll_total,
                      int* ig_nneighbor_all_kind,
                      int* cpu_neighbor,
                      int* ig_cpu_neighbor_info,
                      int* ig_hexa_gnode_glonum, 
                      int* ig_quadp_gnode_glonum,
                      int* ig_quadf_gnode_glonum, 
                      int* ig_hexa_gll_glonum,
                      int* ig_quadp_gll_glonum,
                      int* ig_quadf_gll_glonum, 
                      int* ig_quadp_neighbor_hexa,
                      int* ig_quadp_neighbor_hexaface,
                      int* ig_hexa_material_number)
{ 
   
   mesh_t* mesh = mesh_ptr;
   info_t* info = info_ptr;

   fill_efispec_arrays( info
                           , mesh
                           ,*rank_
                           , ig_ngll_total
                           , ig_nneighbor_all_kind
                           , cpu_neighbor
                           , ig_cpu_neighbor_info
                           , ig_hexa_gnode_glonum
                           , ig_quadp_gnode_glonum
                           , ig_quadf_gnode_glonum
                           , ig_hexa_gll_glonum
                           , ig_quadp_gll_glonum
                           , ig_quadf_gll_glonum
                           , ig_quadp_neighbor_hexa
                           , ig_quadp_neighbor_hexaface
                           , ig_hexa_material_number);

   free_all_2();

   return 0;
}

// write domain files for each proc of efispec3D
void fill_efispec_arrays( info_t* info
                        , mesh_t* mesh
                        , int   rank
                        , int * ig_ngll_total
                        , int * ig_nneighbor_all_kind
                        , int * cpu_neighbor
                        , int * ig_cpu_neighbor_info
                        , int * ig_hexa_gnode_glonum
                        , int * ig_quadp_gnode_glonum
                        , int * ig_quadf_gnode_glonum
                        , int * ig_hexa_gll_glonum
                        , int * ig_quadp_gll_glonum
                        , int * ig_quadf_gll_glonum
                        , int * ig_quadp_neighbor_hexa
                        , int * ig_quadp_neighbor_hexaface
                        , int * ig_hexa_material_number)
{ 
   int iproc, ielt, i, j, k;
   proc_info_t* proc_tab = info->proc;
   elt_info_t*  elt_tab  = info->elt;
      
   int icpu=0;//nb connected proc
   for(iproc = 0; iproc<mesh->npart; iproc++) {
      if (proc_tab->connex[iproc]) cpu_neighbor[icpu++] = iproc; // list of connected proc
   }

   int iparax = 0;
   int ifsurf = 0;

   for (ielt=0; ielt<proc_tab->nb_elt; ielt++) {// for each hexa of proc iproc

      int ihexa = ielt+1;
      int ielt_num, ielt_number, ielt_face, ielt_edge, ielt_corner, ielt_coty, ielt_cpu;

      elt_info_t* curhex = proc_tab->local_elts[ielt];
      int globnum = curhex->globalnum;

      // materiau
      ig_hexa_material_number[ielt] = mesh->layer[globnum];
      // geometry
      for (i=0; i<mesh->num_nodes_hexa; i++) {// geom coords of hexa geom nodes
                                              // use efispec node order
                                              // ig_hexa_gnode_glonum(inode,ihexa)

         if (mesh->eind[mesh->eptr[globnum]+corner2efispec[i]-1] == 0) {
            STOP("num node should nor be 0 (fortran arrays start at 1)");
         }

         ig_hexa_gnode_glonum[ielt * mesh->num_nodes_hexa + i] = mesh->eind[mesh->eptr[globnum]+corner2efispec[i]-1]; // numerotation starts from 1
      }

      mod_init_mesh_mp_init_gll_number_(&ihexa,ig_ngll_total);

      // faces
      for (i=0; i<NFACE; i++) {
         int iface = i+1;
         // 6 faces connectivity
         int type = curhex->faces[i].type;
         switch (type) {
            case HEXA :   // hexa
                     ielt_num = elt_tab[curhex->faces[i].num_neigh].localnum + 1; // num commence à 1
                     ielt_face = face2efispec[curhex->faces[i].num_conn];
                     ielt_coty = curhex->faces[i].orientation;
                     ielt_cpu = mesh->part[curhex->faces[i].num_neigh];
                     if (ielt_cpu == rank) {
                        
                        mod_init_mesh_mp_propagate_gll_nodes_face_(&ihexa, &iface, &ielt_num, &ielt_face, &ielt_coty);
                        
                     } else {
                        if (*ig_nneighbor_all_kind >= 26*info->proc->nb_ext) STOP("ig_nneighbor_all_kind too large\n");

                        ig_cpu_neighbor_info[*ig_nneighbor_all_kind*3 + 0] = ihexa;
                        ig_cpu_neighbor_info[*ig_nneighbor_all_kind*3 + 1] = iface;
                        ig_cpu_neighbor_info[*ig_nneighbor_all_kind*3 + 2] = ielt_cpu;
                        (*ig_nneighbor_all_kind)++;
                     }
                     break;

            case QUAD_P :   // quad p
                     ig_quadp_neighbor_hexa[iparax] = ihexa;
                     ig_quadp_neighbor_hexaface[iparax] = iface;
                     if (mesh->num_node_per_dim*mesh->num_node_per_dim == 4) {   
                        switch(iface) {
                           case 1 : ig_quadp_gnode_glonum[iparax*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+0];
                                    ig_quadp_gnode_glonum[iparax*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+1];
                                    ig_quadp_gnode_glonum[iparax*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+2];
                                    ig_quadp_gnode_glonum[iparax*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+3];
                                    break;
                           case 2 : ig_quadp_gnode_glonum[iparax*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+1];
                                    ig_quadp_gnode_glonum[iparax*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+0];
                                    ig_quadp_gnode_glonum[iparax*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+4];
                                    ig_quadp_gnode_glonum[iparax*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+5];
                                    break;
                           case 3 : ig_quadp_gnode_glonum[iparax*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+2];
                                    ig_quadp_gnode_glonum[iparax*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+1];
                                    ig_quadp_gnode_glonum[iparax*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+5];
                                    ig_quadp_gnode_glonum[iparax*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+6];
                                    break;
                           case 4 : ig_quadp_gnode_glonum[iparax*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+3];
                                    ig_quadp_gnode_glonum[iparax*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+2];
                                    ig_quadp_gnode_glonum[iparax*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+6];
                                    ig_quadp_gnode_glonum[iparax*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+7];
                                    break;
                           case 5 : ig_quadp_gnode_glonum[iparax*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+0];
                                    ig_quadp_gnode_glonum[iparax*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+3];
                                    ig_quadp_gnode_glonum[iparax*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+7];
                                    ig_quadp_gnode_glonum[iparax*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+4];
                                    break;
                           case 6 : ig_quadp_gnode_glonum[iparax*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+7];
                                    ig_quadp_gnode_glonum[iparax*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+6];
                                    ig_quadp_gnode_glonum[iparax*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+5];
                                    ig_quadp_gnode_glonum[iparax*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+4];
                                    break;
                        }
                     } else {
                        STOP("not yet ready for 27 nodes elements\n");
                     }
                     iparax++;
                     
                     mod_init_mesh_mp_propagate_gll_nodes_quad_(&ihexa, &iface, &iparax, ig_quadp_gll_glonum, &info->proc->nb_quad_p);
                     
                     break;
            case QUAD_F :   // quad f
                     if (mesh->num_node_per_dim*mesh->num_node_per_dim == 4) {   
                        switch(iface) {
                           case 1 : ig_quadf_gnode_glonum[ifsurf*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+0];
                                    ig_quadf_gnode_glonum[ifsurf*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+1];
                                    ig_quadf_gnode_glonum[ifsurf*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+2];
                                    ig_quadf_gnode_glonum[ifsurf*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+3];
                                    break;
                           case 2 : ig_quadf_gnode_glonum[ifsurf*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+0];
                                    ig_quadf_gnode_glonum[ifsurf*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+4];
                                    ig_quadf_gnode_glonum[ifsurf*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+5];
                                    ig_quadf_gnode_glonum[ifsurf*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+1];
                                    break;
                           case 3 : ig_quadf_gnode_glonum[ifsurf*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+1];
                                    ig_quadf_gnode_glonum[ifsurf*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+5];
                                    ig_quadf_gnode_glonum[ifsurf*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+6];
                                    ig_quadf_gnode_glonum[ifsurf*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+2];
                                    break;
                           case 4 : ig_quadf_gnode_glonum[ifsurf*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+3];
                                    ig_quadf_gnode_glonum[ifsurf*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+2];
                                    ig_quadf_gnode_glonum[ifsurf*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+6];
                                    ig_quadf_gnode_glonum[ifsurf*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+7];
                                    break;
                           case 5 : ig_quadf_gnode_glonum[ifsurf*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+0];
                                    ig_quadf_gnode_glonum[ifsurf*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+3];
                                    ig_quadf_gnode_glonum[ifsurf*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+7];
                                    ig_quadf_gnode_glonum[ifsurf*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+4];
                                    break;
                           case 6 : ig_quadf_gnode_glonum[ifsurf*4+0] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+4];
                                    ig_quadf_gnode_glonum[ifsurf*4+1] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+7];
                                    ig_quadf_gnode_glonum[ifsurf*4+2] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+6];
                                    ig_quadf_gnode_glonum[ifsurf*4+3] = ig_hexa_gnode_glonum[ielt*mesh->num_nodes_hexa+5];
                                    break;
                        }
                     } else {
                        STOP("not yet ready for 27 nodes elements\n");
                     }
                     ifsurf++;
                     
                     mod_init_mesh_mp_propagate_gll_nodes_quad_(&ihexa, &iface, &ifsurf, ig_quadf_gll_glonum, &info->proc->nb_quad_f);
                     
                     break;
         }
      }
      // edges
      for (i=0; i<NEDGE; i++) {  // 12 edges connectivity
         int iedge = i+1;
         conn_info* ptr = &curhex->edges[i];
         int nb_neighbor = getnbneigh(EDGE, curhex, i);
         for (int ineigh = 0; ineigh < nb_neighbor; ineigh++) {
            switch(ptr->type) {
               case HEXA :    ielt_number = elt_tab[ptr->num_neigh].localnum + 1;   // start from 1
                              ielt_edge   = edge2efispec[ptr->num_conn];            // num edge of neighbor
                              ielt_coty   = ptr->orientation;                       // orientation
                              ielt_cpu    = mesh->part[ptr->num_neigh];             // proc of neighbour
                              if (ielt_cpu == rank) {
                                 if (ielt_number > ihexa) mod_init_mesh_mp_propagate_gll_nodes_edge_(&ihexa, &iedge, &ielt_number, &ielt_edge, &ielt_coty);
                              } else {
                                 ig_cpu_neighbor_info[*ig_nneighbor_all_kind*3 + 0] = ihexa;
                                 ig_cpu_neighbor_info[*ig_nneighbor_all_kind*3 + 1] = NFACE + iedge;
                                 ig_cpu_neighbor_info[*ig_nneighbor_all_kind*3 + 2] = ielt_cpu;
                                 (*ig_nneighbor_all_kind)++;
                              }
               case NONE : break;
               default   : STOP("edge neighbor should be an hexa");
            }
            ptr = ptr->next;
         }
      }
      // corner
      for (j=0; j<NCORNER; j++) { // 8 corners connectivity
         int inode = j+1;
         i = cornerefispec2cubit[j]; // written in efispec order
         conn_info* ptr = &curhex->corners[i];
         int nb_neighbor = getnbneigh(CORNER, curhex, i); // nb of neighbors connected to corner i
         for (k=0; k<nb_neighbor; k++) {
            switch(ptr->type) {
               case HEXA :    ielt_number = elt_tab[ptr->num_neigh].localnum + 1;   // start from 1
                              ielt_corner = corner2efispec[ptr->num_conn];          // num edge of neighbor
                              ielt_cpu    = mesh->part[ptr->num_neigh];             // proc of neighbour
                              if (ielt_cpu == rank) {
                                 if (ielt_number > ihexa) mod_init_mesh_mp_propagate_gll_nodes_corner_(&ihexa, &inode, &ielt_number, &ielt_corner);
                              } else {
                                 ig_cpu_neighbor_info[*ig_nneighbor_all_kind*3 + 0] = ihexa;
                                 ig_cpu_neighbor_info[*ig_nneighbor_all_kind*3 + 1] = NFACE + NEDGE + inode;
                                 ig_cpu_neighbor_info[*ig_nneighbor_all_kind*3 + 2] = ielt_cpu;
                                 (*ig_nneighbor_all_kind)++;
                              }
               case NONE : break;
               default   : STOP("corner neighbor should be an hexa");
            }
            ptr = ptr->next;
         }
      }
   }

   if (iparax != info->proc->nb_quad_p || ifsurf != info->proc->nb_quad_f) STOP("problem with quad number");

   return;
}

// get all topological informations on all elements. return a filled info_t struct
// write geom nodes in files and overwrite mesh->eind with local node num
info_t* getConnInfo(mesh_t* mesh, int rank, int* ig_mesh_nnode)
{ 
   int ielt, iproc, inodeptr;
   
   // data structure allocation
   proc_info_t* proc_tab = MALLOC(proc_info_t, 1);       // DAVID2FREE
   elt_info_t * elt_tab  = MALLOC(elt_info_t, mesh->nh); // DAVID2FREE
   info_t     * info     = MALLOC(info_t, 1);            // DAVID2FREE
   info->proc            = proc_tab;
   info->elt             = elt_tab;
   proc_tab->connex      = MALLOC(char, mesh->npart);    // DAVID2FREE

   // alloc & settings
   memset(proc_tab->connex, 0, mesh->npart*sizeof(char));
   proc_tab->nb_elt    = 0;
   proc_tab->nb_conn   = 0;
   proc_tab->nb_ext    = 0;
   proc_tab->nb_quad_p = 0;
   proc_tab->nb_quad_f = 0;
   
   // get all connectivity info for each element
   for(ielt=0; ielt<mesh->nh; ielt++) {
      //printf("                           \relt : %d / %d", ielt, mesh->nh);
      iproc = mesh->part[ielt];
      if (iproc == rank) {
         elt_tab[ielt].globalnum = ielt;
         
         // check if elem is outer & set proc connectivity info
         elt_tab[ielt].outer = getProcConn(ielt, iproc, proc_tab, mesh);

         if (elt_tab[ielt].outer) proc_tab->nb_ext++;

         proc_tab->nb_elt++;
         // element's neighbors connectivity
         setEltConn(ielt, iproc, elt_tab, mesh, proc_tab);
      }
   }
   
   // count local geom nodes
   idx_t* nodeMask;
   nodeMask = MALLOC(idx_t, mesh->nn);
   memset(nodeMask, 0, mesh->nn*sizeof(idx_t));
   int nb_local_nodes = 0;
   for(ielt=0; ielt<mesh->nh; ielt++) {
      if (mesh->part[ielt] == rank) {
         for(inodeptr = mesh->eptr[ielt]; inodeptr<mesh->eptr[ielt+1]; inodeptr++) {
            if (nodeMask[mesh->eind[inodeptr]] == 0) {
               ++nb_local_nodes;
                                        nodeMask[mesh->eind[inodeptr]]++;
            }
         }   
      }
   }
   *ig_mesh_nnode = nb_local_nodes;

   // allocate arrays
   xcoord = (float*)malloc(nb_local_nodes*sizeof(float));
   ycoord = (float*)malloc(nb_local_nodes*sizeof(float));
   zcoord = (float*)malloc(nb_local_nodes*sizeof(float));
   // fill arrays 
   memset(nodeMask, 0, mesh->nn*sizeof(idx_t));
   nb_local_nodes = 0;
   for(ielt=0; ielt<mesh->nh; ielt++) {
      if (mesh->part[ielt] == rank) {
         for(inodeptr = mesh->eptr[ielt]; inodeptr<mesh->eptr[ielt+1]; inodeptr++) {
            
            int globalnodenum = mesh->eind[inodeptr];
            
            if (nodeMask[globalnodenum] == 0) {
               nodeMask[globalnodenum] = ++nb_local_nodes;
               mesh->eind[inodeptr] = nb_local_nodes;   // overwrite eind since we don't use it anymore

               // warning 0 or 1 C/FORTRAN convention 
               xcoord[nb_local_nodes-1] = mesh->ncoords[globalnodenum].x;
               ycoord[nb_local_nodes-1] = mesh->ncoords[globalnodenum].y;
               zcoord[nb_local_nodes-1] = mesh->ncoords[globalnodenum].z;
            } else {
               int localnumber = nodeMask[globalnodenum];
               mesh->eind[inodeptr] = localnumber;      // overwrite eind since we don't use it anymore
            }
         }   
      }
   }
      assert(*ig_mesh_nnode == nb_local_nodes);
   free(nodeMask);
   
   // alloc & settings
   int num_outer, num_inner;
   num_outer = 0;
   num_inner = proc_tab->nb_ext;
   proc_tab->local_elts = MALLOC(elt_info_t*, proc_tab->nb_elt); // DAVID2FREE
   
   // build data struct for this proc
   for(ielt=0; ielt<mesh->nh; ielt++) {
      int elt_proc = mesh->part[ielt];
      if (elt_proc == rank) {
         if (elt_tab[ielt].outer) {
            proc_tab->local_elts[num_outer] = &elt_tab[ielt];
            elt_tab[ielt].localnum = num_outer++; // warning, localnum start from 0
         } else {
            proc_tab->local_elts[num_inner] = &elt_tab[ielt];
            elt_tab[ielt].localnum = num_inner++; // warning, localnum start from 0
         }
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
      elt_tab[ielt].faces[i].type    = NONE;
      elt_tab[ielt].faces[i].next    = NULL;
      elt_tab[ielt].faces[i].defined = 0;
   }
   for (i=0; i<12; i++) {
      elt_tab[ielt].edges[i].type    = NONE;
      elt_tab[ielt].edges[i].next    = NULL;
      elt_tab[ielt].edges[i].defined = 0;
   }
   for (i=0; i<8; i++) {
      elt_tab[ielt].corners[i].type    = NONE;
      elt_tab[ielt].corners[i].next    = NULL;
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
      if (neigh >= mesh->nh) {   // quad   
         if (mesh->types[neigh] == QUAD_P) {
            proc_tab->nb_quad_p++;
         } else {
            proc_tab->nb_quad_f++;         
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
   int face, inode, inode2;
   
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
         if (!proc_tab->connex[neig_proc]) {
            proc_tab->connex[neig_proc] = 1;
            proc_tab->nb_conn++;
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


void readCubitMesh_1(char* fileinp, mesh_t* mesh, int* number_of_elemnt_block)
{ 
   char line[LENMAX];
   int iblock = 0;
   FILE *unit_inp;
   char *hexa8=NULL, *hexa27=NULL;
   int readquad = 0;

   unit_inp = fopen(fileinp,"r");
   if(unit_inp == NULL)
   {
     STOP("can't open input file (readCubitMesh_1)")
   }

   //First pass on file *.inp: count number of nodes, hexa, quad_parax and quad_fsurf
   mesh->nh = mesh->nq_parax = mesh->nq_surf = 0;
   while (!feof(unit_inp))
   {
      char * cdum = fgets(line, LENMAX, unit_inp);

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
   fclose(unit_inp);
   return;
}
 
// read a cubit mesh and fill mesh struct
void readCubitMesh_2(char* fileinp, mesh_t* mesh, int* number_of_elemnt_block)
{ 
   char  line[LENMAX];
   int   ielt, inode, itmp;
   int   loc_eind = 0;
   int   ielt_glo = 0;
   int   iblock   = 0;
   FILE *unit_inp;
   char *hexa8=NULL, *hexa27=NULL;
   char *quadf=NULL, *quadp=NULL;
   int   ihexa=0;

   unit_inp = fopen(fileinp,"r");
   if(unit_inp == NULL)
   {
     STOP("can't open input file (readCubitMesh_2)")
   }
  
    //Second pass on file *.inp: make eind & eptr
   iblock = ihexa = 0;
   while (!feof(unit_inp))
   {   
      char *cdum = fgets(line, LENMAX, unit_inp);

      //fill array gnodes with coordinates of geometric nodes
      if (strstr(line,"NSET=ALLNODES"))
         for (inode=0;inode<mesh->nn;inode++) {
            int idum = fscanf(unit_inp,"%d," "%"SCREAL "," "%"SCREAL "," "%"SCREAL, &itmp, &mesh->ncoords[inode].x, &mesh->ncoords[inode].y, &mesh->ncoords[inode].z); //SCREAL is defined in metis.h
         }

      //fill eptr and eind for hexa
      if ((hexa8 = strstr(line,HEXA_8_HDR)) || (hexa27 = strstr(line,HEXA_27_HDR)))
      {
         // get num of layer
         char *hdr = hexa8?HEXA_8_HDR:HEXA_27_HDR;
         int layer = atoi(rtrim(strstr(line,hdr) + strlen(hdr)));
         
         for (ielt=0;ielt<number_of_elemnt_block[iblock];ielt++)
         {
            int idum = fscanf(unit_inp,"%d",&itmp);
            mesh->eptr[ielt_glo]    = loc_eind;
            mesh->types[ielt_glo++] = HEXA;
            mesh->layer[ihexa++]    = layer;
            for (inode=0;inode<mesh->num_nodes_hexa;inode++)
            {
               int idum = fscanf(unit_inp,",%d",&mesh->eind[loc_eind]);
               mesh->eind[loc_eind++]--;   // because num starts on 1 in cubit meshes
            }
         }
         
         iblock++;
      }
      
       //fill eptr and eind for quad_parax and quad_fsurf
      if ((quadp = strstr(line,QUAD_P_HDR)) || (quadf = strstr(line,QUAD_F_HDR)))
      {
         for (ielt=0;ielt<number_of_elemnt_block[iblock];ielt++)
         {
            int idum = fscanf(unit_inp,"%d",&itmp);
            mesh->eptr[ielt_glo]    = loc_eind;
            mesh->types[ielt_glo++] = quadp?QUAD_P:QUAD_F;
            for (inode=0;inode<N_NODES_QUAD;inode++)
            {
               int idum =  fscanf(unit_inp,",%d",&mesh->eind[loc_eind]);
               mesh->eind[loc_eind++]--;   // because num starts on 1 in cubit meshes
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
   // sign :   + : left->right or bottom->up
   //         - : left<-right or bottom<-up
   // ie ortn[] = {1, -3} ->
   //   +-----> i
   //   |
   //   | k
   //   V
   
   char Hsign, Vsign;

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

   
   if (Vrot&1) ret_or[1] = Vsign*abs(V);   // impair -> Vert
   else ret_or[0] = Vsign*abs(V);         // pair -> Hor

   if (Hrot&1) ret_or[1] = Hsign*abs(H);   // impair -> Vert
   else ret_or[0] = Hsign*abs(H);         // pair -> Hor

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
   
   return 0;   // opposite directions
}


void add_conn(topo_t typecon, elt_info_t* elt_tab, int hexa_src, int numcon_src, idx_t num_neigh, int numcon_neigh, int orientation, elt_t type)
{ 
   conn_info *ptr, *last;

   switch(typecon) {
      case FACE :   if (elt_tab[hexa_src].faces[numcon_src].defined) { // malloc
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
      case FACE :   printf("WARNING : funct getnbneigh() should not be called for faces\n");
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
