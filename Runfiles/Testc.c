# include <stdio.h>
# include <stdlib.h>
# include <cs.h>
# include "clayer.h"
# include <time.h> 

/*
Just a small test routine for compilation etc 
*/
int main(void) 
{
    /* Construct a matrix */
    srand(time(NULL));
    cs *T, *B, *C, *thisCS, *A, *thismmCS; 
    MyCSparse *this;
    MyCSparse thismm;
    int const nval = 500;
    int rows[nval];
    int cols[nval];
    double vals[nval]; 
    int i;
    int ncol;
    int nz; 
    int p1, p2, length, p; 
    int temp;


    for(i = 0; i < nval; i++)
    {
        rows[i] = rand() % 1000;
        cols[i] = rand() % 1000;
        vals[i] = rand() % 1000;
    }

    /* Allocate */
    this = ConstructMyCSparse(1000, 1000, nval, rows, cols, vals);
    //this.nrow = 10;
    //this.ncol = 10;
    //this.nval = 5;
    //this.row = rows;
    //this.col = cols;
    //this.val = vals;

    T = cs_spalloc(1000, 1000, nval, 1, 1);
    B = cs_spalloc(this->nrow, this->ncol, this->nval, 1, 1);
    thismm = SpMM(this, this);

    thismmCS = ConvertMyCSparseToCS(&thismm);
    
    cs_print(thismmCS, i);



    /* Add entries */
    for (i = 0; i < nval; i++)
    {
        cs_entry(T, rows[i], cols[i], vals[i]);
    };

    /* Print */
    //cs_print(T, i);

    /* Try some multiplications*/
    A = cs_compress(T); // convert to csc format for multiplication
    B = cs_multiply(A, A);

    /* Convert to triplet format again - fully based on 
    umfpack_col_to_triplet
     */
    C = B;
    ncol = B->n;
    C->nz = B->p[ncol];
    for (i = 0; i < ncol; i++)
    {
        p1 = B->p[i];
        p2 = B->p[i+1];
        length = p2-p1;

        for(p = p1; p < p2; p++)
        {
            C->p[p] = i;
        }

    }
    //cs_print(C, i);


    return 0;
};


