// // Minimal illustrative MPI implementation of Totoni et al. (2012) structure‑adaptive sparse
// // triangular solve.  The goal is to show how to reproduce the three heuristics without Charm++.
// // * Column‑block distribution (round‑robin).
// // * Analysis step:   – mark independent rows  – reorder diagonal block  – inspect off‑diag rows.
// // * Solve step:      – eager non‑blocking sends to hide communication.
// //
// // Caveats
// // =======
// // 1.  This is a teaching / demo program (<1000 sloc).  It is **not** fully tuned, nor does it handle
// //     every corner case (e.g. zero pivots, very small matrices).
// // 2.  Only single right‑hand side (vector b) is supported.
// // 3.  Input format = simple ASCII CSR identical to the format used在你与 ChatGPT 的对话里:  
// //        n
// //        row_ptr.size()  (== n+1)
// //        row_ptr[0] … row_ptr[n]
// //        col_idx.size()
// //        col_idx[…]
// //        val.size()
// //        val[…]
// //     All indices assumed **0‑based**.
// // 4.  Each rank owns one column‑block (easy to extend to over‑decomposition by splitting further).
// //
// // Author:  Adapted from Totoni, Heath & Kale, SC’12 (structure‑adaptive SpTRSV)
// // ---------------------------------------------------------------------------
// #include <mpi.h>
// #include <vector>
// #include <fstream>
// #include <iostream>
// #include <cassert>
// #include <algorithm>
// #include <numeric>
// #include <chrono>
// #include <sstream>
// #include <iomanip>

// #include "../include/common.h"

// // ---------------------------------------------------------------------------
// CSRMatrix loadCSR(std::string& filename) {
//     std::ifstream infile(filename);
//     if (!infile.is_open()) {
//         throw std::runtime_error("Cannot open CSR file: " + filename);
//     }
//     auto start = std::chrono::high_resolution_clock::now();

//     CSRMatrix csr;
//     int row_ptr_size, col_id_size, val_size;

//     infile >> csr.n;
//     infile >> row_ptr_size;
//     csr.row_ptr.resize(row_ptr_size);
//     for (int i = 0; i < row_ptr_size; ++i) {
//         infile >> csr.row_ptr[i];
//     }

//     infile >> col_id_size;
//     csr.col_id.resize(col_id_size);
//     for (int i = 0; i < col_id_size; ++i) {
//         infile >> csr.col_id[i];
//     }

//     infile >> val_size;
//     csr.val.resize(val_size);
//     for (int i = 0; i < val_size; ++i) {
//         infile >> csr.val[i];
//     }
//     infile.close();

//     auto end = std::chrono::high_resolution_clock::now();
//     double duration = std::chrono::duration<double>(end - start).count();
//     std::cout << "CSR matrix loaded from \"" << filename << "\" in " 
//               << duration << " seconds." << std::endl;
              
//     return csr;
// }
// // ---------------------------------------------------------------------------
// // ---------------------------------------------------------------------------
// inline int ownerOf(int idx,int blkSize){return idx/blkSize;}

// struct LocalBlk
// {
//     int nGlob, blkSize;
//     int colBeg, colEnd;
//     // diagonal rows (within [colBeg,colEnd))
//     std::vector<int> rows;           // global row id per local row
//     std::vector<int> rp, ci;         // CSR (global col indices)
//     std::vector<double> a;
//     // dependency info
//     std::vector<std::vector<int>> dependsOn, dependents;
//     std::vector<int> order;          // execution sequence
//     std::vector<double> diag;        // diagonal entries (same size as rows)
// };

// // ===========  build local block & dense‑offdiag scan  ========================
// static constexpr double TILE_THR = 0.30;   // >30% off‑diag density triggers tiling
// static constexpr int    BUF_CAP  = 64;     // aggregation size

// LocalBlk buildLocal(const CSRMatrix &G,int rank,int P)
// {
//     LocalBlk B;  B.nGlob=G.n;  B.blkSize=(G.n+P-1)/P;
//     B.colBeg = rank*B.blkSize;
//     B.colEnd = std::min(G.n,B.colBeg+B.blkSize);

//     std::vector<int> mapRow2Local(G.n,-1);
//     // pass1: collect rows with any nnz in my column range
//     for(int i=0;i<G.n;++i){
//         bool touched=false;
//         for(int k=G.row_ptr[i];k<G.row_ptr[i+1];++k){
//             int j=G.col_id[k];
//             if(j>=B.colBeg && j<B.colEnd){touched=true;break;}
//         }
//         if(touched){ mapRow2Local[i]=(int)B.rows.size(); B.rows.push_back(i);}    }

//     int m=B.rows.size(); B.rp.assign(m+1,0);

//     size_t offDiagNNZ=0;
//     for(int lr=0;lr<m;++lr){
//         int i=B.rows[lr];
//         for(int k=G.row_ptr[i];k<G.row_ptr[i+1];++k){
//             B.ci.push_back(G.col_id[k]);
//             B.a .push_back(G.val[k]);
//             if(G.col_id[k]>=B.colEnd) ++offDiagNNZ;              // 统计右侧nnz
//         }
//         B.rp[lr+1]=(int)B.ci.size();
//     }

//     // decide dense‑offdiag strategy
//     double dens = double(offDiagNNZ)/double(B.ci.size());
//     if(dens > TILE_THR){
//         if(rank==0) std::cerr<<"[Info] rank"<<rank<<" triggers tile mode (density="<<dens<<")\n";
//         // *demo version* — we do **not** create extra chares; we simply keep the nnz but mark
//         // them so that they get computed as soon as x[i] available.  Real tile splitting needs
//         // extra blocks; omitted for brevity.
//     }

//     B.dependsOn.resize(m);
//     B.dependents.resize(m);
//     B.diag.assign(m,0.0);
//     return B;
// }

// // ===========  analysis: independents & dep graph  ===========================
// void analyze(LocalBlk &B,int rank,int P)
// {
//     int m=B.rows.size();
//     std::vector<char> dep(m,0);
//     for(int lr=0;lr<m;++lr){
//         int i=B.rows[lr];
//         for(int k=B.rp[lr];k<B.rp[lr+1];++k){
//             int j=B.ci[k]; if(j==i) B.diag[lr]=B.a[k];
//             if(j<i){                    // dependency edge i←j
//                 int owner=ownerOf(j,B.blkSize);
//                 if(owner!=rank) dep[lr]=1;              // cross‑block dep
//                 else {
//                     // intra block; will be discovered once dep[j] known later
//                 }
//             }
//         }
//     }
//     // second pass mark intra‑block dep
//     for(int lr=0;lr<m;++lr){
//         int i=B.rows[lr];
//         for(int k=B.rp[lr];k<B.rp[lr+1];++k){
//             int j=B.ci[k]; if(j<i && ownerOf(j,B.blkSize)==rank){
//                 int jr = std::lower_bound(B.rows.begin(),B.rows.end(),j)-B.rows.begin();
//                 if(jr<m && B.rows[jr]==j && dep[jr]) dep[lr]=1;
//             }
//         }
//     }
//     // build execution order (indep rows reversed first)
//     for(int lr=m-1; lr>=0; --lr) if(!dep[lr]) B.order.push_back(lr);
//     for(int lr=0 ; lr< m; ++lr) if( dep[lr]) B.order.push_back(lr);
// }

// // ===========  aggregated send helper ========================================
// struct SendBuf{
//     std::vector<double> payload; std::vector<int> tag; // tag=row id
// };

// void flush(int dst, SendBuf &buf, MPI_Comm comm){
//     if(buf.payload.empty()) return;
//     int n = (int)buf.payload.size();
//     MPI_Send(&n,1,MPI_INT,dst,0,comm);
//     MPI_Send(buf.tag.data(),n,MPI_INT,dst,1,comm);
//     MPI_Send(buf.payload.data(),n,MPI_DOUBLE,dst,2,comm);
//     buf.payload.clear(); buf.tag.clear();
// }

// // ===========  solve kernel  ==================================================
// void sptrsv(LocalBlk &B, std::vector<double>& x, MPI_Comm comm)
// {
//     int rank,P; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&P);
//     int m=B.rows.size();
//     x.assign(B.nGlob,0.0);

//     // aggregated buffers per dst
//     std::vector<SendBuf> out(P);

//     auto sendVal=[&](int dst,int row,double val){
//         out[dst].payload.push_back(val);
//         out[dst].tag    .push_back(row);
//         if((int)out[dst].payload.size()>=BUF_CAP) flush(dst,out[dst],comm);
//     };

//     // === recursion over execution order ===
//     for(int idx=0; idx<m; ++idx){
//         int lr = B.order[idx]; int i=B.rows[lr];

//         // probe for any incoming bulk first (message priority)
//         int flag=1; MPI_Status st;
//         while(flag){
//             MPI_Iprobe(MPI_ANY_SOURCE,0,comm,&flag,&st);
//             if(!flag) break;
//             int cnt; MPI_Recv(&cnt,1,MPI_INT,st.MPI_SOURCE,0,comm,&st);
//             std::vector<int> tag(cnt); std::vector<double> val(cnt);
//             MPI_Recv(tag.data(),cnt,MPI_INT,st.MPI_SOURCE,1,comm,&st);
//             MPI_Recv(val.data(),cnt,MPI_DOUBLE,st.MPI_SOURCE,2,comm,&st);
//             for(int t=0;t<cnt;++t) x[tag[t]]=val[t];
//         }

//         // wait for required remote predecessors (blocking‑probe loop)
//         for(int k=B.rp[lr]; k<B.rp[lr+1]; ++k){
//             int j=B.ci[k]; if(j<i && ownerOf(j,B.blkSize)!=rank){
//                 while(x[j]==0.0){      // crude spin w/ probe
//                     MPI_Iprobe(MPI_ANY_SOURCE,0,comm,&flag,&st);
//                     if(flag){
//                         int cnt; MPI_Recv(&cnt,1,MPI_INT,st.MPI_SOURCE,0,comm,&st);
//                         std::vector<int> tag(cnt); std::vector<double> val(cnt);
//                         MPI_Recv(tag.data(),cnt,MPI_INT,st.MPI_SOURCE,1,comm,&st);
//                         MPI_Recv(val.data(),cnt,MPI_DOUBLE,st.MPI_SOURCE,2,comm,&st);
//                         for(int t=0;t<cnt;++t) x[tag[t]]=val[t];
//                     }
//                 }
//             }
//         }

//         // compute row i
//         double sum=0.0;
//         for(int k=B.rp[lr]; k<B.rp[lr+1]; ++k){ int j=B.ci[k]; if(j<i) sum+=B.a[k]*x[j]; }
//         x[i]=(1.0-sum)/B.diag[lr];   // RHS=1

//         // eager aggregated send to all procs that own rows > i with nnz in col i
//         for(int p=0;p<P;++p){ if(p==rank) continue; sendVal(p,i,x[i]); }
//     }

//     // flush remaining buffers
//     for(int p=0;p<P;++p) flush(p,out[p],comm);
// }

// // ====================  Utils  ===============================================
// template<class V> void bcast_vec(std::vector<V>&v,int root,MPI_Comm comm){
//     int sz=(int)v.size(); MPI_Bcast(&sz,1,MPI_INT,root,comm);
//     if((int)v.size()!=sz) v.resize(sz);
//     MPI_Datatype dt = std::is_same<V,int>::value?MPI_INT:MPI_DOUBLE;
//     MPI_Bcast(v.data(),sz,dt,root,comm);
// }

// // ====================  main  ================================================
// int main(int argc,char**argv){
//     MPI_Init(&argc,&argv);
//     int rank,P; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&P);

//     if(argc<2){ if(rank==0) std::cerr<<"Usage: ./sptrsv matrix.csr\n"; MPI_Finalize(); return 0; }

//     CSRMatrix A;
//     if (rank == 0) 
//     {
//         std::string CSR_file = argv[1];
//         std::ifstream infile(CSR_file);
//         if (!infile) {
//             std::cerr << "Cannot open file: " << CSR_file << std::endl;
//             MPI_Abort(MPI_COMM_WORLD, -1);
//         }
//         std::cout << "Loading A in CSR format..." << std::endl;
//         A = loadCSR(CSR_file);
//         std::cout << "[rank0]  Matrix loaded: n = "
//                   << A.n << ", nnz = " << A.val.size() << '\n';
//     }
//     MPI_Bcast(&A.n,1,MPI_INT,0,MPI_COMM_WORLD);
//     bcast_vec(A.row_ptr,0,MPI_COMM_WORLD);
//     bcast_vec(A.col_id ,0,MPI_COMM_WORLD);
//     bcast_vec(A.val    ,0,MPI_COMM_WORLD);

//     LocalBlk B = buildLocal(A,rank,P);
//     analyze(B,rank,P);

//     std::vector<double> x;
//     double t0=MPI_Wtime();
//     sptrsv(B,x,MPI_COMM_WORLD);
//     double t1=MPI_Wtime();

//     if(rank==0) std::cout<<"Elapsed="<<std::fixed<<std::setprecision(6)<<t1-t0<<" s\n";
//     MPI_Finalize();
//     return 0;
// }
