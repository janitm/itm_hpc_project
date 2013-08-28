#include "common.h"
#include "functions.h"
#include "minunit.h"

static int nn = 0;  //number of cells in one dimension

void test_setup() {
n = 100;
}

void test_teardown() {
// Nothing
}

/*those are just examples from previous project
MU_TEST(test_insertion_sort) {
    //check if the order of the elements in the sorted array is correct
    int i, condition;
    float_t data[n];

    //create unsorted array and set condition to 0
    for(i=n; i > 0; i--)
        data[i]= 0.1 * (float)i;
    condition = 0;

    //sort data and check for correct order
    insertion_sort(data,n);
    for(i = 0; i < (n-1); i++)
    {
        if( (data[i+1] - data[i]) < 0. )
            condition++;
    }
    mu_assert(condition == 0, "insertion sort is not working correctly: wrong order");

    //recreate unsorted array and set condition to 0 again
    for(i=n; i > 0; i--) data[i]= 0.1 * (float)i;
    condition = 0;

    //check for the right number of elements in data
    condition = (int) sizeof(data) / sizeof(data[0]);
    mu_assert(condition == n, "insertion sort is not working correctly: wrong number of elements");

}

MU_TEST(test_fill_matrix)
{
    int i,j, condition;
    float_t **matrix;
    float_t data[n];

    for(i=0; i<n; i++)
        data[i] = 1./((float) i+1);

    matrix = calloc(n, sizeof(float_t **));
    for(i = 0; i < n; i++)
        matrix[i] = calloc(n, sizeof(float_t *));

    fill_matrix(data,matrix,n);

    //go through matrix and check correct fill up
    condition = 0;
    for(i = 0; i <n-1; i++)
        for(j = 0; j < n-1; j++)
            if(matrix[i][j+1] > matrix[i][j] || matrix[i+1][j] > matrix[i][j])
                condition += 1;
    mu_assert(condition == 0, "matrix fill up did not work correctly");
    //printf("last matrix el: %f \n", matrix[n-1][n-1]);
    //mu_assert(matrix[n-1][n-1] == 1./(float_t)((n+1)*(n+1)), "matrix fill up did not work correctly");
}


MU_TEST(test_generate_histogram)
{
    int i,j, boxes, bin_sum, check_sum;
    float_t **matrix;
    float_t max, check_max;
    int* histogram;

    bin_sum = 0;
    max = 0.;
    matrix = calloc(n, sizeof(float_t **));
    boxes = n / 2;
    histogram = malloc(boxes * sizeof(int));


    //allocate matrix
    for(i = 0; i < n; i++)
       matrix[i] = calloc(n, sizeof(float_t *));

    //generate test matrix
    for(i = 0; i < n; i++)
    {
        for(j=0; j<n; j++)
        matrix[i][j] = 1./n * ((float)i) * 1./n * ((float)j); //0.1*((float_t)i) * 0.1*((float_t)j);
    }//end for

    //test the histogram function
    max = generate_histogram(matrix, histogram, n, boxes);

    //sum over all entries in histogram
    for(i = 0; i<boxes; i++)
        bin_sum = bin_sum + histogram[i];

    //sum of all entries in histogram should equal n*n -2*n - 2*(n-2)
    check_sum = (n*n) - 2*n - 2*(n-2);
    printf("\n check_sum: %d", check_sum);
    printf("\n bin_sum: %d \n", bin_sum);
    mu_assert(bin_sum == ((n*n) - 2*n - 2*(n-2)), "wrong number of elements in histogram");
    free(matrix);
    free(histogram);

    /*
    //check maximum element, since in this construction it is given by expression below
    //C arrays start with 0, therefore substract 1
    free(matrix);
     //allocate matrix
    for(i = 0; i < n; i++)
       matrix[i] = calloc(n, sizeof(float_t *));

    //generate test matrix
    for(i = 0; i < n; i++)
        for(j=0; j<n; j++)
            matrix[i][j] = (float_t) i*j;
    max = generate_histogram_v2(matrix, histogram, n, boxes);

    n -= 1;
    check_max = (float_t)( fabsf((float)((n-1)*(n-1) - (n-1)*(n)))
                            +fabsf((float)((n-1)*(n-1) - (n-1)*(n-2)))
                            +fabsf((float)((n-1)*(n-1) - (n)*(n-1)))
                            +fabsf((float)((n-1)*(n-1) - (n-2)*(n-1))));
    //cast to integer in assert, since both values differ slighlty due to numerics
    check_max /= (float_t) n;
    printf("\ncheck_max: %f\n", check_max);
    printf("\n max: %f \n", max);
    printf("\n diff: %f \n", fabsf(check_max - max));

    mu_assert( fabsf(check_max - max) <= 0.1, "wrong maximum element in histogram");*/
}
*/


MU_TEST_SUITE(test_suite) {
MU_SUITE_CONFIGURE(&test_setup, &test_teardown);

MU_RUN_TEST(test_insertion_sort);
MU_RUN_TEST(test_fill_matrix);
MU_RUN_TEST(test_generate_histogram);

}

int main(int argc, char *argv[]) {
MU_RUN_SUITE(test_suite);
MU_REPORT();
return 0;
}
