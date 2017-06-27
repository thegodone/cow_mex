#include <stdlib.h>
#include "mex.h"
#include <vector>
#include <cmath>


struct IndexCoeff
{
    double *Coeff;
    int *Index;
};
struct IndPosValue
{
    double* inds;
    double* posValue;
    int len;
};

bool debug=false;

int getMaxPosition(double* X, int n){
    int maxPos = 0;
    for (int i=1; i< n; i++){
        if (X[i]> X[maxPos]) maxPos=i;
    }
    return maxPos;
}
double* getCoeff_b(double* DA, int n,int* Allowed_Arcs_index, int N_AA){
    double* Coeff_b = (double*) mxMalloc (n*N_AA * sizeof(double)); //new double[n*N_AA];
    int col_i;
    for (int k=0; k < N_AA; k++){
        col_i = Allowed_Arcs_index[k];
        for (int i=0; i<n; i++){
            Coeff_b[k*n +i] = DA[col_i*n +i]  ;//be changed to DA[col_i*n + i]
        }
    }
    return Coeff_b;
}

double* getElementMatStyle(double* DA,  int* index, int index_length, int element_len){
    
    double* data = (double*) mxMalloc (index_length * sizeof(double)); //new double[index_length];
    for (int i=0; i<index_length; i++){
        int el_at = index[i]-1;
        if (el_at < element_len)
            data[i] = DA[el_at] ;
        else
            data[i] = 0;
    }
    return data;
}
int* getElementMatStyle(int* DA,  int* index, int index_length, int element_len){
    
    int* data = (int*) mxMalloc (index_length * sizeof(int)); // new int[index_length];
    for (int i=0; i<index_length; i++){
        int el_at = index[i]-1;
        if (el_at <element_len)
            data[i] = DA[el_at] ;
        else
            data[i]=0;
    }
    return data;
}

double* SubArray(double* A, int* B, int n){
    double* res = (double*) mxMalloc (n * sizeof(double));// new double[n];
    
    for (int i=0; i<n; i++){
        res[i] = A[i] - B[i];
    }
    return res;
}
int* addToMatrix(int* data, int value, int n){
    int* res = (int*) mxMalloc (n * sizeof(int));//new int[n];
    
    for (int i=0; i<n; i++){
        res[i] = data[i] + value;
    }
    return res;
}

int* getIndex_Node(int* DA, int valueSub,  int n, int* Allowed_Arcs_index, int N_AA, int value){
    int* index_Nodes = (int*) mxMalloc (n*N_AA * sizeof(int));//new int[n*N_AA];
    int col_i;
    for (int k=0; k < N_AA; k++){
        col_i = Allowed_Arcs_index[k];
        for (int i=0; i<n; i++){
            // DA[col_i*n +i]
            index_Nodes[k*n +i] = DA[col_i*n + i] + value -1 -valueSub;// ajust value according to C++
              
        }
    }
    
    return index_Nodes;
}

double* modifyXiSeg(double* Xi_Seg, double*Xi_diff, double* Coeff_b, int  n ){
    
    double* res = (double*) mxMalloc (n * sizeof(double)); //new double[n];
    for (int i=0; i<n; i++){
        res[i] = Xi_Seg[i] + Coeff_b[i] * Xi_diff[i];
    }
    return res;
}
double* getXi_Seg_mean(double* Xi_Seg,  int  n, int N_AA ){
    double*  Xi_Seg_mean = (double*) mxMalloc (N_AA * sizeof(double)); // new double[N_AA];
    for (int k=0; k < N_AA; k++){
        Xi_Seg_mean[k] =0;
        for (int i=0; i<n; i++){
            Xi_Seg_mean[k] += Xi_Seg[k*n + i];
        }
        Xi_Seg_mean[k] = Xi_Seg_mean[k] / n;
    }
    return Xi_Seg_mean;
}

double* getNorm_Xi_Seg_cen(double* Xi_Seg, double* Xi_Seg_mean, int  n, int N_AA ){
    double Xi_sum_sqr ;
    double* Norm_Xi_Seg_cen =(double*) mxMalloc (N_AA * sizeof(double)); // new double[N_AA];
    for (int k=0; k < N_AA; k++){
        Xi_sum_sqr =0;
        for (int i=0; i<n; i++){
            Xi_sum_sqr += Xi_Seg[k*n + i]*Xi_Seg[k*n + i] ;
        }
        
        Norm_Xi_Seg_cen[k] = sqrt(std::abs(Xi_sum_sqr  - n * Xi_Seg_mean[k]*Xi_Seg_mean[k])); // need L2 cannot remove sqrt!
    }
    
    return    Norm_Xi_Seg_cen;
}
double* getCCs_Values(double* Xi_Seg, double *TSeg_centred,int  n, int N_AA){
    double* CCs_Values = (double*) mxMalloc (N_AA * sizeof(double));// new double[N_AA];
    for (int k=0; k < N_AA; k++){
        CCs_Values[k] = 0;
        
        for (int i=0; i<n; i++){ 
            CCs_Values[k] += TSeg_centred[i] * Xi_Seg[k*n + i];
        }
    }
    
    return CCs_Values;
}
double* getFactor(double* Xi_Seg_mean, double* Norm_Xi_Seg_cen,  double Norm_TSeg_cen, int N_AA){
    double* factor = (double*) mxMalloc (N_AA * sizeof(double)); //new double[N_AA];
    for (int k=0; k < N_AA; k++){
        factor[k] = Norm_TSeg_cen * Norm_Xi_Seg_cen[k];
    }
    return factor;
}
double* getCOS_FUN(double* table_2_values , double* CCs_Values, double* factors, int N_AA){
    
    double* Cost_Fun = (double*) mxMalloc (N_AA * sizeof(double)); //new double[N_AA];
    for (int k=0; k < N_AA; k++){
        double table_val = table_2_values[k];
        double factor = factors[k];
        double CCs_Node = 0;
        if (factor > 0){
            
            CCs_Node = CCs_Values[k] / factor;
        }
        Cost_Fun[k] = CCs_Node + table_val;
    }
    
    return Cost_Fun;
    
}



int* calculatePrec_Node(int table_val, int* a, int n){
    int* Prec_Node = (int*) mxMalloc (n * sizeof(int)); //new int[n];
    for (int i=0; i<n; i++){
        
        Prec_Node[i] = table_val - (int)a[i];
    }
    return Prec_Node;
}
int* calculateAllowed_Arcs_index(int* Allowed_Arcs, int Prec_Nodes_length,  int N_AA){
    int* Allowed_Arcs_index = (int*) mxMalloc (N_AA * sizeof(int)); //new int[N_AA];
    int j=0;
    for (int i=0; i<Prec_Nodes_length; i++){
        if (Allowed_Arcs[i]==1){
            Allowed_Arcs_index[j] = i;
            j++;
        }
    }
    return Allowed_Arcs_index;
}
int* calculateNodes_TablePointer(int* Allowed_Arcs,int* Prec_Nodes, int b, int Prec_Nodes_length,  int N_AA){
    int* Nodes_TablePointer = (int*) mxMalloc (N_AA * sizeof(int)); //new int[N_AA];
    
    int j=0;
    for (int i=0; i<Prec_Nodes_length; i++){
        if (Allowed_Arcs[i]==1){
            Nodes_TablePointer[j]= b + Prec_Nodes[i];
            j++;
        }
    }
    return Nodes_TablePointer;
}

int calculateAllowed_N_AA(int* Allowed_Arcs, int Prec_Nodes_length){
    int N_AA=0;
    for (int i=0; i<Prec_Nodes_length; i++){
        if (Allowed_Arcs[i]==1){
            N_AA++;
        }
    }
    return N_AA;
}
int* calculateAllowed_ArcsNodes(int* Prec_Nodes, int Prec_Nodes_length, int boundry1,int boundry2){
    int * Allowed_Arcs = (int*) mxMalloc (Prec_Nodes_length * sizeof(int));// new int[Prec_Nodes_length];
    for (int i=0; i<Prec_Nodes_length; i++){
        if (Prec_Nodes[i]>=boundry1 && Prec_Nodes[i]<=boundry2){
            Allowed_Arcs[i]=1;
        }else{
            Allowed_Arcs[i]=0;
        }
    }
    
    return Allowed_Arcs;
}
double norm(double *X, int n){
    
    double m_sum = 0;
    for (int i=0; i<n;++i)
    {
        m_sum  += X[i]* X[i];
    }
    
    return sqrt(m_sum);
}

double* calculateTg_center(double* T, int start, int len){
    
    double* res = (double*) mxMalloc (len * sizeof(double)); //new double[len];
    double sum =0;
    for (int i=0; i<len; i++){
        res[i] = T[i+start];
        sum += res[i] /len;
    }
    for (int i=0; i<len; i++){
        res[i] = res[i] -sum;
    }
    return res;
    
}

mxArray * getMexDoubleArray(double* v, int v_len){
    mxArray * mx =  mxCreateNumericMatrix(1,v_len, mxDOUBLE_CLASS, mxREAL);
    double *indexes = mxGetPr(mx);
    for (int i=0; i<v_len; i++){
        indexes[i] = v[i];
    }
    mxFree(v);
    return mx;
}



int* addToIntArray(int* X, int val, int n){
    int* IntX = (int*) mxMalloc (n * sizeof(int)); //new int[n];
    for (int i=0;i<n;i++)
        IntX[i] =  X[i] + val;
    return IntX;
}
int* getLenSeg(int nSeg,int seg, int T_len, int X_len){
    int* LengSeg = (int*) mxMalloc (2*nSeg * sizeof(int)); //new int[2*nSeg];
   
    for (int i=0; i<nSeg; i++){
        
        LengSeg[i] = seg -1; //LenSeg[i] 
        LengSeg[nSeg + i] = seg -1; //LenSeg[i]
    }
    LengSeg[nSeg-1] += (T_len-1) % LengSeg[0];//LenSeg[nSeg-1]
    LengSeg[2*nSeg-1] += (X_len-1) % LengSeg[nSeg];//LenSeg[nSeg-1]
    return LengSeg;
}
int* getSlacks_vec(int slack, int n){
    int *V = (int*) mxMalloc (n * sizeof(int)); //new int[n];
    int initial = -slack;
    for (int i=0; i< n; i++){
        V[i] = initial + i;
    }
    return V;
    
}
int* getBTORBP(int* LengSegRow, int start, int nSeg){
    int* V = (int*) mxMalloc ((nSeg+1) * sizeof(int)); //new int[nSeg+1];
    V[0] =1;
    for (int i=1; i<=nSeg;i++)
        V[i] = V[i-1] + LengSegRow[start+i-1];
    return V;
}

int* getBounds(int* bP, int slack, int nSeg){
    int* bounds = (int*) mxMalloc ((2*nSeg+2) * sizeof(int)); // new int[2*nSeg+2];
   
    for (int i=0; i<=nSeg;i++){
        int bound_a= bP[i] +i*-1*slack;
        int bound_a_2= bP[i] + i*slack;
        int bound_b = bP[i] + (nSeg-i) *-1*slack; //i=0 , -824
        int bound_b_2 = bP[i] + (nSeg-i)*slack;//i=0 , 826
		
        bounds[i] = bound_a; //i=0 , 1 //bounds[i]
        bounds[nSeg + 1 +i] =bound_a_2; //i=0 , 1 //
        
        
        if (bounds[i]<bound_b)bounds[i] =bound_b; //bounds[i]
        if (bounds[nSeg + 1 +i]>bound_b_2)bounds[nSeg + 1 +i] =bound_b_2;//bounds[nSeg+1 +i]
    }
    
    return bounds;
}
double* diff(double* X, int nX){
    double* Xdiff = (double*) mxMalloc ((nX-1) * sizeof(double)); //new double[nX-1];
    for (int i=1; i<nX;i++){
        Xdiff[i-1] = X[i] - X[i-1] ;
    }
    return Xdiff;
}

int* getTableIndex(int* bounds, int bound_len ){ //int bounds
    int* Table_Index = (int*) mxMalloc ((bound_len+1) * sizeof(int)); //new int [bound_len+1];
    
    Table_Index[0] = 0;
    for (int i=1; i<=bound_len;i++)
        Table_Index[i] = Table_Index[i-1 ]  + bounds[bound_len +i-1 ] - bounds[i-1] +1;
    
    return Table_Index;
}

double* getTable(int* bounds, int table_len,int bounds_Len){
    
    double* Table = (double*) mxMalloc ((3*table_len) * sizeof(double)); //new double[3*table_len];//never use third row
   
    int i;
    
    Table[table_len +0] =0;
    Table[2*table_len +0] =0;
    for (i=1;i<table_len;i++){
        Table[table_len + i] =-1e10; /* -infinity */
        Table[2*table_len + i] =0;
    }
    int k=0;
    int val;
    //mexPrintf("Table_Index:");
    for (i=0; i<bounds_Len; i++){
        for (val=bounds[i];val<=bounds[bounds_Len +i];val++){ //from bounds[i] bounds[bounds_Len + i]
            Table[k] = val;
            k++;
            //mexPrintf("%f ", round(val));
            //if (k%25 == 0) mexPrintf("\n ");
        }
    }
    //mexPrintf("\n ");
    return Table;
}
int* getArrayFromTo(int from, int to){
    int n= to - from +1;
    int* X = (int*) mxMalloc (n * sizeof(int)); //new int[n];
    for (int i=0;i<n;i++ ){
        X[i] = from +i;
    }
    return X;
}

double* getPArray(double q, int t){
    int n = q +1;
    double* P = (double*) mxMalloc (n * sizeof(double)); // new double[n];
    for (int i=0;i<n;i++ ){
        P[i] = ((double)(i * t ))/q + 1;
    }
    return P;
}
int* caclulateK(double* P, int P_len, int* PP, int PP_len){
    
    int* K = (int*) mxMalloc (P_len * sizeof(int)); //new int[P_len];
    int j=0;
    for (int i=0; i<P_len; i++) {
        
        while ( j< PP_len ) {
            
            if (P[i]>=PP[j]){
                
                j++;
            }
            else{
                
                break ;
            }
            
        }
        K[i] =j;
        //mexPrintf("i %d, j %d  PP_len %d K[i] %d \n", i, j,PP_len, K[i]);
        if (K[i]<1) K[i]=1;
        if (K[i]>=PP_len)  K[i] = PP_len - 1;
        //mexPrintf("i %d, j %d  PP_len %d K[i] %d \n", i, j,PP_len, K[i]);
    }
    return K;
}

IndexCoeff getIndexCoeff(int n, int LenSeg_bound, int* Slacks_vec, int Slacks_vec_len){
    int* nPrime = addToMatrix(Slacks_vec,LenSeg_bound , Slacks_vec_len);
 
    int A_m = Slacks_vec_len;
    int A_len =  n; //LenSeg[0][0]+1;
    double* A = (double*) mxMalloc ((Slacks_vec_len * n) * sizeof(double)); //new double[Slacks_vec_len * n];
    int B_m =  Slacks_vec_len;
    int B_len =n; //LenSeg[0][0]+1;
    
    int* B = (int*) mxMalloc ((Slacks_vec_len * n) * sizeof(int)); //new int[Slacks_vec_len * n];
    
    int* PP, PP_len;
    double* P;
    int P_len;
    int* K;
    P_len = A_len;
    for (int i=0; i<Slacks_vec_len; i++){
        //mexPrintf("loop i %d :\n", i);
        PP = getArrayFromTo(1, nPrime[i]);
        PP_len = nPrime[i];
        //mexPrintf("PP:\n");
        //printArray(PP, PP_len);
        P = getPArray(A_len-1, nPrime[i]-1);
        //mexPrintf("P:\n");
        //printDoubleArray(P, P_len);
        K= caclulateK(P, P_len, PP, PP_len);
        
        //mexPrintf("K:\n");
        //printArray(K, P_len);
        
        int* pp_k =getElementMatStyle(PP,  K, P_len, PP_len );
        
        //mexPrintf("pp K:\n");
        //mexPrintf("i=%d\n", i);
        //printArray(pp_k, P_len);
        for (int k=0; k<P_len; k++){
			A[i*P_len +k] = P[k] - pp_k[k];
            B[i*P_len +k] = K[k] -Slacks_vec[i];
		}
      
        mxFree(K),
        mxFree(P);
        mxFree(PP);
        mxFree(pp_k);
    }
    
    mxFree(nPrime);
    IndexCoeff AB;
    AB.Coeff = A;
    AB.Index = B;
    return AB;
}
IndPosValue caculateCow(double *T,int T_len ,  double *X, int X_len, int seg, int slack ){
    //accept seg
    int nSeg =  floor((T_len - 1) / (seg - 1));
//    mexPrintf("nSeg = %d\n", nSeg);
    
    
    int Slacks_vec_len = 2*slack+1;
    int *Slacks_vec = getSlacks_vec(slack, Slacks_vec_len);
    
    //printArray(Slacks_vec, Slacks_vec_len);
    
    
    int LenSeg_Len = nSeg; //mxGetN(prhs[2]);
    int LenSeg_m =   2; //mxGetM(prhs[2]);
    int* LenSeg= getLenSeg( nSeg, seg,  T_len,  X_len); //getIntTwoDimArray(mxGetPr(prhs[2]), LenSeg_Len, LenSeg_m);
    
    
    int* bT = getBTORBP(LenSeg,0, nSeg);
    int bT_len = nSeg+1;
    //printArray(bT, bT_len);
    int* bP = getBTORBP(LenSeg, nSeg, nSeg);
    //mexPrintf("bP: \n  ");
    ////printArray(bP, bT_len);
    int* Bounds = getBounds(bP,  slack,  nSeg);
    int Bounds_Len = nSeg +1 ;
    int Bounds_m =  2;
    double *Xdiff = diff( X, X_len);//array
    int Xdiff_Len = X_len-1;
    int Table_Index_len = Bounds_Len +1;

    int *Table_Index = getTableIndex( Bounds,  Bounds_Len );
    
    
    int Table_Len = Table_Index[Table_Index_len-1];
    
    double* Table =getTable(Bounds, Table_Len,Bounds_Len);
    //mexPrintf("Table first row: len %d, first %f, last %f\n", Table_Len,Table[0], Table[Table_Len-1] );
    
    IndexCoeff AB = getIndexCoeff(LenSeg[0]+1, LenSeg[nSeg]+1, Slacks_vec, Slacks_vec_len);//LenSeg[0]+1, LenSeg[LenSeg_Len]+1
    
    //mexPrintf("AB ");
    int A_m = Slacks_vec_len;
    
    
    int A_len =  LenSeg[0]+1;
    double *A = AB.Coeff;
    int B_m =  Slacks_vec_len;
    int B_len = LenSeg[0]+1;
    
    int  *B = AB.Index;
    
    AB = getIndexCoeff( LenSeg[nSeg-1]+1, LenSeg[2*nSeg-1]+1, Slacks_vec, Slacks_vec_len); //LenSeg[nSeg-1]+1, LenSeg[2*nSeg-1]+1
    int A_last_len = LenSeg[nSeg-1]+1;
    int A_last_m =  Slacks_vec_len;
    double *A_last = AB.Coeff;
    
    
    int B_last_len = LenSeg[nSeg-1]+1;
    int B_last_m =  Slacks_vec_len;
    int *B_last =  AB.Index;
    
    
    int N_AA, pos_value;
    double ind;
    
    
    for (int i_seg = 0; i_seg<nSeg; i_seg++){
        
        //mexPrintf("i_seg %d\n", i_seg);
        int* a = addToIntArray(Slacks_vec, LenSeg[nSeg + i_seg], Slacks_vec_len);
        //printArray(a,Slacks_vec_len);
        
        
        int b             = Table_Index[i_seg]+ 1 - Bounds[i_seg];
        int c             = LenSeg[i_seg] + 1;
        int count  = 0;
        int Node_Z        = Table_Index[i_seg + 2];                     // Last node for segment i_seg
        int Node_A        = Table_Index[i_seg+1];                // First node for segment i_seg
        
        int Int_Index_seg_len ,Int_Coeff_Seg_len;
        int *Int_Index_Seg;
        double *Int_Coeff_Seg;
        if (i_seg < (nSeg-1)){
            Int_Index_Seg = B; //addToMatrix(B, (0- LenSeg[1][i_seg]),  A_m,A_len); //modifcation needed
            Int_Coeff_Seg = A;
            Int_Index_seg_len =B_len ;
            Int_Coeff_Seg_len = A_len;
            
        }else {
            Int_Index_Seg = B_last; //addToMatrix(B_last , (0- LenSeg[1][i_seg]), A_last_m, A_last_len); //modification needed
            Int_Coeff_Seg = A_last;
            Int_Index_seg_len = B_last_len ;
            Int_Coeff_Seg_len = A_last_len;
        }
        
        int TSeg_centred_len = bT[i_seg + 1] - bT[i_seg] +1;
        //mexPrintf("TSeg_centred_len %d\n", TSeg_centred_len);
        double *TSeg_centred =  calculateTg_center(T, bT[i_seg]-1, TSeg_centred_len);
        //printDoubleArray(TSeg_centred, TSeg_centred_len);
        double Norm_TSeg_cen = norm(TSeg_centred, TSeg_centred_len);
        for (int i_node= Node_A; i_node< Node_Z; i_node++){
            // if(i_node == 3188) debug=true;
            //     else
            debug= false;
            double tab_val = Table[i_node];
            
            int Prec_Nodes_length = Slacks_vec_len;
            int  boundry1 = Bounds[i_seg];
            int boundry2 = Bounds[Bounds_Len +i_seg];
            //mexPrintf("boundry1 %d, boundry2 %d\n", boundry1, boundry2);
            int  ind = 0;
            double pos_value = 0;
            
            
            int* Prec_Nodes = calculatePrec_Node(tab_val,  a,Prec_Nodes_length);
            
            int* Allowed_Arcs  = calculateAllowed_ArcsNodes(Prec_Nodes,  Prec_Nodes_length,  boundry1, boundry2);
            //mexPrintf("Prec_Nodes");
            //printArray(Prec_Nodes,Prec_Nodes_length );
            int N_AA = calculateAllowed_N_AA( Allowed_Arcs,  Prec_Nodes_length);
            //mexPrintf("i_node %d, N_AA %d\n", i_node, N_AA);
            Table[Table_Len + i_node] = 0;
            Table[2*Table_Len + i_node] = 0;
            
            if (N_AA>0){
                
                int Nodes_TablePointer_len = N_AA;
                int* Allowed_Arcs_index =  calculateAllowed_Arcs_index(Allowed_Arcs,  Prec_Nodes_length,   N_AA);
                
                int* Nodes_TablePointer = calculateNodes_TablePointer( Allowed_Arcs,Prec_Nodes, b,  Prec_Nodes_length,   N_AA);
                
                int* Index_Node = getIndex_Node(Int_Index_Seg, 
                                                LenSeg[nSeg +i_seg], Int_Index_seg_len, Allowed_Arcs_index, N_AA, tab_val);//need c
                
                double* Coeff_b = getCoeff_b(Int_Coeff_Seg,  Int_Coeff_Seg_len, Allowed_Arcs_index,  N_AA); //need change 
                
                double* Xi_diff = getElementMatStyle(Xdiff,Index_Node, Int_Index_seg_len*N_AA, Xdiff_Len );
                double* Xi_Seg_tmp  =   getElementMatStyle(X,Index_Node, Int_Index_seg_len*N_AA, X_len );
                
                double* Xi_Seg = modifyXiSeg(Xi_Seg_tmp, Xi_diff, Coeff_b,Int_Index_seg_len*N_AA);
                
                double* Xi_Seg_mean = getXi_Seg_mean( Xi_Seg,  Int_Index_seg_len, N_AA );
                double*  Norm_Xi_Seg_cen = getNorm_Xi_Seg_cen( Xi_Seg,  Xi_Seg_mean, Int_Index_seg_len, N_AA);
                
                double* table_2_values = (double*) mxMalloc (N_AA * sizeof(double)); //new double[N_AA];
                for (int k=0; k < N_AA; k++){
                    int index_Nodes_TablePointer = (int)Nodes_TablePointer[k];
                    table_2_values[k]= Table[Table_Len +index_Nodes_TablePointer-1];
                }
                
                double* CCs_values = getCCs_Values(Xi_Seg, TSeg_centred,  Int_Index_seg_len,  N_AA);
                
                double* factors = getFactor( Xi_Seg_mean, Norm_Xi_Seg_cen, Norm_TSeg_cen, N_AA);
                
                double* Cost_Fun = getCOS_FUN(table_2_values , CCs_values, factors, N_AA);
                
                
                int pos =  getMaxPosition(Cost_Fun, N_AA);
                
                Table[ Table_Len + i_node] =  Cost_Fun[pos];
                Table[2*Table_Len + i_node] =  Nodes_TablePointer[pos];
                mxFree(Allowed_Arcs_index);
                mxFree(Nodes_TablePointer);
                mxFree(Index_Node);
                mxFree(Coeff_b);
                mxFree(Xi_diff);
                mxFree(Xi_Seg);
                mxFree(Xi_Seg_tmp);
                mxFree(Xi_Seg_mean);
                mxFree(Norm_Xi_Seg_cen);
                mxFree(table_2_values);
                mxFree(CCs_values);
                mxFree(factors);
                mxFree(Cost_Fun);
            }
            //mexPrintf("modifying table");
            
            mxFree(Prec_Nodes);
            mxFree(Allowed_Arcs);
            count++;
        }
        
       // mxFree(Int_Index_Seg;//should be modified
        mxFree(a);
        mxFree(TSeg_centred);
        
    }
    
    IndPosValue results ;
    results.inds = (double*) mxMalloc (Table_Len * sizeof(double)); //new double[Table_Len]; ////;
    results.posValue = (double*) mxMalloc (Table_Len * sizeof(double)); //new double[Table_Len];
    results.len = Table_Len;
    
    for (int i=0; i<Table_Len; i++){
        results.inds[i] = Table[Table_Len + i];
        results.posValue[i] = Table[2*Table_Len + i];
    }
    
     mxFree(Slacks_vec);
     mxFree(LenSeg);//
     mxFree(bP);
     mxFree(bT);
     mxFree(Bounds);//
     mxFree(Xdiff);
     mxFree(Table_Index);
     mxFree(Table);//
     mxFree(A);
     mxFree(B);
     mxFree(A_last);
     mxFree(B_last);
    return results;
    
}

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
    
    double *T = mxGetPr(prhs[0]);
    int T_len = mxGetN(prhs[0]);
    
    double *X = mxGetPr(prhs[1]); //array
    int X_len = mxGetN(prhs[1]);
    int seg = mxGetScalar(prhs[2]);
    int slack = mxGetScalar(prhs[3]);
    
    IndPosValue results = caculateCow(T,T_len ,  X,  X_len,  seg,  slack );
    
    
    
    plhs[0] = getMexDoubleArray(results.inds, results.len);
    
    plhs[1] = getMexDoubleArray(results.posValue, results.len);
    
    
}
