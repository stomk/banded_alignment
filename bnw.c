#include <stdlib.h>
#include <stdio.h>

typedef struct{
    char* seq1_aln;
    char* seq2_aln;
    int aln_size;
    int score;
} alignment;


int max2(int a, int b){
    return a >= b ? a : b;
}

int max3(int a, int b, int c){
    int max = a;
    if(b > max) max = b;
    if(c > max) max = c;
    return max;
}

void free_alignment(alignment* aln){
    if(aln->seq1_aln){ free(aln->seq1_aln); aln->seq1_aln = NULL;}
    if(aln->seq2_aln){ free(aln->seq2_aln); aln->seq2_aln = NULL;}
    if(aln){ free(aln); aln = NULL; }
}

void print_matrix(int** A, int n, int m, int k){
    int d = m - n;
    int i,j;
    for(i = 0; i <= k; i++){
        printf("%3d   ", i);
        for(j = 0; j <= i+k+d+1; j++)
            printf("%d  ", A[i][j]);
        printf("\n");
    }
    for(i = k+1; i <= n-k-1; i++){
        printf("%3d   ", i);
        for(j = 0; j <= 2*k+d+2; j++)
            printf("%d  ", A[i][j]);
        printf("\n");
    }
    for(i = n-k; i <= n; i++){
        printf("%3d   ", i);
        for(j = 0; j <= -i+k+d+n+1; j++)
            printf("%d  ", A[i][j]);
        printf("\n");
    }
}


alignment* align(char* seq1, int seq1_len, char* seq2, int seq2_len,
                 int band_tolerance,
                 int match_score, int mismatch_score,
                 int gap_open_penalty, int gap_extend_penalty){
    int n;
    int m;
    char* a = NULL;
    char* b = NULL;
    int** M = NULL;
    int** I = NULL;
    int** J = NULL;
    int** H = NULL;
    int k = band_tolerance / 2;
    int d;
    int e = match_score;
    int f = mismatch_score;
    int g = gap_open_penalty;
    int h = gap_extend_penalty;
    int i, j;
    int l;
    char* x = NULL;
    char* y = NULL;
    char* p = NULL;
    char* q = NULL;
    int max_aln_size;
    alignment* aln = NULL;
    int aln_size;

    
    // Set sequences for alignment
    if(seq1_len < seq2_len){
        n = seq1_len;
        m = seq2_len;
        a = calloc(n + 2, sizeof(char));
        b = calloc(m + 2, sizeof(char));
        for(i = 1; i <= n; i++){
            a[i] = seq1[i-1];
        }
        for(j = 1; j <= m; j++){
            b[j] = seq2[j-1];
        }
    }else{
        n = seq2_len;
        m = seq1_len;
        a = calloc(n + 2, sizeof(char));
        b = calloc(m + 2, sizeof(char));
        for(i = 1; i <= n; i++){
            a[i] = seq2[i-1];
        }
        for(i = 1; i <= m; i++){
            b[i] = seq1[i-1];
        }
    }

    // distance between lengths of two sequences
    d = m - n;

    // Allocate memory for matrices
    M = (int**)malloc(sizeof(int*) * (n + 1));
    I = (int**)malloc(sizeof(int*) * (n + 1));
    J = (int**)malloc(sizeof(int*) * (n + 1));
    H = (int**)malloc(sizeof(int*) * (n + 1));

    for(i = 0; i <= k; i++){
        M[i] = (int*)malloc(sizeof(int) * (i + k + d + 2));
        I[i] = (int*)malloc(sizeof(int) * (i + k + d + 2));
        J[i] = (int*)malloc(sizeof(int) * (i + k + d + 2));
        H[i] = (int*)malloc(sizeof(int) * (i + k + d + 2));
    }
    for(i = k+1; i <= n-k-1; i++){
        M[i] = (int*)malloc(sizeof(int) * (2*k + d + 3));
        I[i] = (int*)malloc(sizeof(int) * (2*k + d + 3));
        J[i] = (int*)malloc(sizeof(int) * (2*k + d + 3));
        H[i] = (int*)malloc(sizeof(int) * (2*k + d + 3));
    }
    for(i = n-k; i <= n; i++){
        M[i] = (int*)malloc(sizeof(int) * (-i + k + d + n + 2));
        I[i] = (int*)malloc(sizeof(int) * (-i + k + d + n + 2));
        J[i] = (int*)malloc(sizeof(int) * (-i + k + d + n + 2));
        H[i] = (int*)malloc(sizeof(int) * (-i + k + d + n + 2));
    }


    /* Initialize matrices */
    const int z = -10000000;

    // i = 0
    M[0][0] = 0;
    I[0][0] = z;
    J[0][0] = z;
    H[0][0] = 0;

    for(j = 1; j <= k+d+1; j++){
        M[0][j] = z;
        I[0][j] = z;
        J[0][j] = g + h*(j-1);
        H[0][j] = J[0][j];
    }

    // 1 <= i <= k
    for(i = 1; i <= k; i++){
        M[i][0] = z;
        M[i][i + k + d + 1] = z;
        I[i][0] = g + h*(i-1);
        I[i][i + k + d + 1] = z;
        J[i][0] = z;
        H[i][0] = I[i][0];
    }

    // k+1 <= i <= n-k-1
    for(i = k+1; i <= n-k-1; i++){
        M[i][0] = z;
        M[i][2*k + d + 2] = z;
        I[i][2*k + d + 2] = z;
        J[i][0] = z;
    }

    // n-k <= i <= n
    for(i = n-k; i <= n; i++){
        M[i][0] = z;
        J[i][0] = z;
    }

    /* End of Initiallization */
    

    /* Recursion */
    // 1 <= i <= k+1
    for(i = 1; i <= k+1; i++){
        for(j = 1; j <= k+d+i; j++){
            M[i][j] = max3( M[i-1][j-1], I[i-1][j-1], J[i-1][j-1] )
                + (a[i] == b[j] ? e : f);
            I[i][j] = max2( M[i-1][j] + g, I[i-1][j] + h) ;
            J[i][j] = max2( M[i][j-1] + g, J[i][j-1] + h );
            H[i][j] = max3( M[i][j], I[i][j], J[i][j] );
        }
    }

    // k+2 <= i <= n-k
    for(i = k+2; i <= n-k; i++){
        for(j = 1; j <= 2*k+d+1; j++){
            M[i][j] = max3(M[i-1][j], I[i-1][j], J[i-1][j] )
                + (a[i] == b[j+i-k-1] ? e : f);
            I[i][j] = max2( M[i-1][j+1] + g, I[i-1][j+1] + h );
            J[i][j] = max2( M[i][j-1] + g, J[i][j-1] + h );
            H[i][j] = max3( M[i][j], I[i][j], J[i][j] );
        }
    }

    // n-k+1 <= i <= n
    for(i = n-k+1; i <= n; i++){
        for(j = 1; j <= n+k+d+1-i; j++){
            M[i][j] = max3(M[i-1][j], I[i-1][j], J[i-1][j] )
                + (a[i] == b[j+i-k-1] ? e : f);
            I[i][j] = max2( M[i-1][j+1] + g, I[i-1][j+1] + h );
            J[i][j] = max2( M[i][j-1] + g, J[i][j-1] + h );
            H[i][j] = max3( M[i][j], I[i][j], J[i][j] );
        }
    }

    /* End of Recursion */

    
    /* Traceback and fill alignment strings */
    //    max_aln_size = 2*m + 2*k - n + 10;
    max_aln_size = 2*m;

    x = (char*)calloc(max_aln_size + 1, sizeof(char));
    y = (char*)calloc(max_aln_size + 1, sizeof(char));
    p = &x[max_aln_size - 1];
    q = &y[max_aln_size - 1];

    i = n; j = k+d+1;
    int score = H[i][j];
    while(i >= k+2){
        l = j + i - k -1;
        if(H[i][j] == I[i][j]){
            *p = a[i];
            *q = '-';
            i--; j++;
        }else if(H[i][j] == J[i][j]){
            *p = '-';
            *q = b[l];
            j--;
        }else{
            *p = a[i];
            *q = b[l];
            i--;
        }
        p--; q--;
    }
    while( i > 0 || j > 0 ){
        if( H[i][j] == I[i][j] ){
            *p = a[i];
            *q = '-';
            i--;
        }else if( H[i][j] == J[i][j] ){
            *p = '-';
            *q = b[j];
            j--;
        }else{
            *p = a[i];
            *q = b[j];
            i--; j--;
        }
        p--; q--;
    }

    /* End of Traceback */

    // Shift chars in x and y to the head
    i = 0;
    while( x[i] == 0 ) i++;
    aln_size = max_aln_size - i;
    if(i != 0){
        for(p = &x[0], q = &y[0] ; i < max_aln_size; i++){
            *(p++) = x[i];
            *(q++) = y[i];
        }        
        *p = '\0';
        *q = '\0';
    }
        

    // Free memory
    if(a){ free(a); a = NULL;}
    if(b){ free(b); b = NULL;}
    for(i = 0; i <= n; i++){
        if(M[i]){ free(M[i]); M[i] = NULL;}
        if(I[i]){ free(I[i]); I[i] = NULL;}
        if(J[i]){ free(J[i]); J[i] = NULL;}
        if(H[i]){ free(H[i]); H[i] = NULL;}
    }
    if(M){ free(M); M = NULL;}
    if(I){ free(I); I = NULL;}
    if(J){ free(J); J = NULL;}
    if(H){ free(H); H = NULL;}


    // Return alignment
    aln = (alignment*)calloc(1, sizeof(alignment*));
    if(seq1_len < seq2_len){
        aln->seq1_aln = x;
        aln->seq2_aln = y;
    }else{
        aln->seq1_aln = y;
        aln->seq2_aln = x;
    }
    aln->aln_size = aln_size;
    aln->score = score;

    return aln;
}


/* int main(int argc, char**argv){ */
/*     char s1[172] = "GAATTCTCAGTAACTTCTTTGTGTTTGTGTGTATTCAACTCAACAGAGTTGAACTTTCTTTAGAGAGACGCAGAGTTGAAAACCTCTGTTTTGGAATTTGCAAGTGCAGATTTCAAGCGCTTCCTAGCCTATGGCAGAAAAGGAAATATCTTCGTATAAAAACTACACAGA"; */
/*     char s2[172] = "TCATTCTCAACAACTACTTTGTGATGTGTGCGTTCAACTCACAGAGTTTAACCTTTCTTTTCGTAGAGCAGTTTGGAAACACTCTGTTTGTAAAGCCTGGCAAGTGCTTTTTTGGACTTCATTGAGCGCTTCGTTGGAGAAACGGGATTTCTTCATATAATGCTAGACAGA"; */
    
/*     alignment* aln = align(s1, 171, s2, 171, 10, 1, -2, -4, -1); */
/*     printf("%s\n%s\n%d\n", aln->seq1_aln, aln->seq2_aln, aln->aln_size); */
/*     free_alignment(aln); */
    
/*     return 0; */
/* } */
