/*-------------------------------------------------------------
 * build_L_2level.c
 * Create a 160×160 red-black grid, build the 5-point Laplacian,
 * permute (red first, black second), and keep the (≤-diag) part
 * as CSR.  Output: three files row_ptr.txt, col_idx.txt, val.txt
 *------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h> // -lm

#define MAX_NZ_PER_ROW 5         /* 5-point stencil          */

typedef struct {                 /* simple COO triplet       */
    int32_t r, c;
    double  v;
} nz_t;

/* map 2-D (ix,iy) → global id */
static inline int gid(int ix,int iy, int NX)
{ 
    return iy*NX + ix; 
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return EXIT_FAILURE;
    }

    int N = atoi(argv[1]);
    int NX = (int)(floor(sqrt((double)N) + 0.5));
    if (NX * NX != N) {
        fprintf(stderr, "Error: n = %d is not a perfect square\n", N);
        return EXIT_FAILURE;
    }
    int NY = NX;
    int NN = NX * NY;  /* = N */

    /* 2. 在 COOmatrix 中生成 5 点 Laplacian */
    size_t max_nz = (size_t)NN * 5;
    nz_t *coo = malloc(max_nz * sizeof(nz_t));
    if (!coo) { perror("malloc coo"); return EXIT_FAILURE; }
    size_t nz = 0;
    for (int iy = 0; iy < NY; iy++) {
        for (int ix = 0; ix < NX; ix++) {
            int p = gid(ix, iy, NX);
            /* 自身 */
            coo[nz++] = (nz_t){p, p, 4.0};
            /* 四邻域 */
            if (ix > 0)     coo[nz++] = (nz_t){p, gid(ix-1, iy, NX), -1.0};
            if (ix < NX-1)  coo[nz++] = (nz_t){p, gid(ix+1, iy, NX), -1.0};
            if (iy > 0)     coo[nz++] = (nz_t){p, gid(ix, iy-1, NX), -1.0};
            if (iy < NY-1)  coo[nz++] = (nz_t){p, gid(ix, iy+1, NX), -1.0};
        }
    }

    /* 3. 红黑重排 —— 得到 perm 和 invp */
    int32_t *perm = malloc(NN * sizeof(int32_t));
    int32_t *invp = malloc(NN * sizeof(int32_t));
    if (!perm || !invp) { perror("malloc perm"); return EXIT_FAILURE; }
    /* ------ 统计黑点、红点个数 ---------------------------- */
    int black_cnt = 0;
    for (int iy = 0; iy < NY; ++iy)
        for (int ix = 0; ix < NX; ++ix)
            if ((ix + iy) & 1) ++black_cnt;

    int red_cnt   = NN - black_cnt;   /* 红 = 总 - 黑    */
    int rptr = 0;                     /* 红块从 0 开始   */
    int bptr = red_cnt;               /* 黑块紧随其后   */

    /* ------ 建 permutation --------------------------------- */
    for (int iy = 0; iy < NY; ++iy)
        for (int ix = 0; ix < NX; ++ix) 
        {
            int p = gid(ix, iy, NX);
            if (((ix + iy) & 1) == 0) {         /* red */
                perm[p] = rptr;
                invp[rptr++] = p;
            } else {                            /* black */
                perm[p] = bptr;
                invp[bptr++] = p;
            }
        }


    /* 4. 统计并组装下三角 CSR */
    int32_t *row_ptr = calloc(NN+1, sizeof(int32_t));
    if (!row_ptr) { perror("calloc row_ptr"); return EXIT_FAILURE; }
    /* 第一遍计数 */
    for (size_t k = 0; k < nz; k++) {
        int i = perm[coo[k].r], j = perm[coo[k].c];
        if (j <= i) row_ptr[i+1]++;
    }
    for (int i = 0; i < NN; i++)
        row_ptr[i+1] += row_ptr[i];
    size_t nzL = row_ptr[NN];

    int32_t *col_id = malloc(nzL * sizeof(int32_t));
    double  *val    = malloc(nzL * sizeof(double));
    if (!col_id || !val) { perror("malloc col_id/val"); return EXIT_FAILURE; }
    int32_t *cursor = malloc(NN * sizeof(int32_t));
    for (int i = 0; i < NN; i++) cursor[i] = row_ptr[i];
    /* 第二遍填充 */
    for (size_t k = 0; k < nz; k++) {
        int i = perm[coo[k].r], j = perm[coo[k].c];
        double v = coo[k].v;
        if (j <= i) {
            int pos = cursor[i]++;
            col_id[pos] = j;
            val[pos]    = v;
        }
    }

    /* Output format:
        * n
        * row_ptr.size()
        * row_ptr
        * col_id.size()
        * col_id
        * val.size()
        * val
    */
   char fname[64];
   snprintf(fname, sizeof(fname), "syn_matrix_n%d.txt", N);
   FILE *fp = fopen(fname, "w");
   if (!fp) { perror("fopen output file"); return EXIT_FAILURE; }


   fprintf(fp, "%d\n", NN);
   fprintf(fp, "%d\n", NN+1);
   for (int i = 0; i <= NN; i++)
       fprintf(fp, "%d%c", row_ptr[i], (i==NN?'\n':' '));
   fprintf(fp, "%zu\n", nzL);
   for (size_t k = 0; k < nzL; k++)
       fprintf(fp, "%d%c", col_id[k], (k+1==nzL?'\n':' '));
   fprintf(fp, "%zu\n", nzL);
   for (size_t k = 0; k < nzL; k++)
       fprintf(fp, "%.17g%c", val[k], (k+1==nzL?'\n':' '));

   fclose(fp);
   printf("Finished: %s （%d×%d in Lower CSR，%zu nnz）\n",
          fname, NN, NN, nzL);

   free(coo);
   free(perm);
   free(invp);
   free(row_ptr);
   free(col_id);
   free(val);
   free(cursor);

   return EXIT_SUCCESS;
}
