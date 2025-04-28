/*-------------------------------------------------------------
 * build_L_2level.c
 * Create a 160×160 red-black grid, build the 5-point Laplacian,
 * permute (red first, black second), and keep the (≤-diag) part
 * as CSR.  Output: three files row_ptr.txt, col_idx.txt, val.txt
 *------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define NX 160
#define NY 160
#define NN (NX*NY)               /* 25 600 */

#define MAX_NZ_PER_ROW 5         /* 5-point stencil          */

typedef struct {                 /* simple COO triplet       */
    int32_t r, c;
    double  v;
} nz_t;

/* map 2-D (ix,iy) → global id */
static inline int gid(int ix,int iy){ return iy*NX + ix; }


int main(void)
{
    /* 1. build Laplacian in COO ------------------------------*/
    nz_t *coo = malloc(NN*MAX_NZ_PER_ROW*sizeof(nz_t));
    size_t nz = 0;

    for (int iy=0; iy<NY; ++iy)
    for (int ix=0; ix<NX; ++ix){
        int p = gid(ix,iy);

        /* centre */
        coo[nz++] = (nz_t){p,p, 4.0};

        /* four neighbours if inside domain */
        if(ix>0   ) coo[nz++] = (nz_t){p,gid(ix-1,iy),-1.0};
        if(ix<NX-1) coo[nz++] = (nz_t){p,gid(ix+1,iy),-1.0};
        if(iy>0   ) coo[nz++] = (nz_t){p,gid(ix,iy-1),-1.0};
        if(iy<NY-1) coo[nz++] = (nz_t){p,gid(ix,iy+1),-1.0};
    }

    /* 2. create red-black permutation ------------------------*/
    int32_t *perm = malloc(NN*sizeof(int32_t));
    int32_t *invp = malloc(NN*sizeof(int32_t));
    int rptr = 0, bptr = 0;
    /* first count reds to know where black block starts */
    for(int iy=0;iy<NY;++iy)
    for(int ix=0;ix<NX;++ix)
        if( (ix+iy)&1 ) ++bptr;
    rptr = 0;                /* reds at front         */
    int bbase = bptr;        /* black block offset    */

    for(int iy=0;iy<NY;++iy)
    for(int ix=0;ix<NX;++ix){
        int p = gid(ix,iy);
        if( ((ix+iy)&1)==0 ){            /* red */
            perm[p] = rptr;
            invp[rptr++] = p;
        }else{                           /* black */
            perm[p] = bptr;
            invp[bptr++] = p;
        }
    }

    /* 3. convert to CSR keeping only col<=row (lower part) ----*/
    int32_t *row_ptr = calloc(NN+1, sizeof(int32_t));
    int32_t *col_idx = malloc(nz * sizeof(int32_t));
    double  *val     = malloc(nz * sizeof(double));
    size_t   nzL = 0;

    /* first pass: count non-zeros per row of L                */
    for(size_t k=0;k<nz;++k){
        int32_t i = perm[coo[k].r];
        int32_t j = perm[coo[k].c];
        if(j<=i){ ++row_ptr[i+1]; }
    }
    for(int i=0;i<NN;i++) row_ptr[i+1]+=row_ptr[i];

    /* second pass: fill arrays                                */
    int32_t *cursor = malloc(NN*sizeof(int32_t));
    for(int i=0;i<NN;i++) cursor[i]=row_ptr[i];

    for(size_t k=0;k<nz;++k){
        int32_t i = perm[coo[k].r];
        int32_t j = perm[coo[k].c];
        if(j<=i){
            size_t pos = cursor[i]++;
            col_idx[pos] = j;
            val[pos]     = coo[k].v;
            ++nzL;
        }
    }

    /* 4. write to matrix.txt —— */
    /* Output format:
        * n
        * row_ptr.size()
        * row_ptr
        * col_id.size()
        * col_id
        * val.size()
        * val
        */
    
    FILE *fp = fopen("synthetic_matrix.txt","w");
  
    /* n */
    fprintf(fp, "%d\n", NN);

    /* row_ptr.size() */
    fprintf(fp, "%d\n", NN+1);
    /* row_ptr */
    for (int i = 0; i <= NN; i++) {
        fprintf(fp, "%d%s", row_ptr[i], (i==NN ? "\n" : " "));
    }

    /* col_id.size() */
    fprintf(fp, "%zu\n", nzL);
    /* col_id */
    for (size_t k = 0; k < nzL; k++) {
        fprintf(fp, "%d%s", col_idx[k], (k+1==nzL ? "\n" : " "));
    }

    /* val.size() */
    fprintf(fp, "%zu\n", nzL);
    /* val */
    for (size_t k = 0; k < nzL; k++) {
        fprintf(fp, "%.17g%s", val[k], (k+1==nzL ? "\n" : " "));
    }

    fclose(fp);
    printf("CSR matrix file saved.\n");

    /* tidy */
    free(coo); free(perm); free(invp);
    free(row_ptr); free(col_idx); free(val); free(cursor);

    
    return 0;
}
