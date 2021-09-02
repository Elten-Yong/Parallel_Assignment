#include <iostream>
#include <omp.h>
#include <mpi.h>
#include <assert.h>
#define NUM_THREADS 8
#define lastNumber 300000000
using namespace std;


int main2()
{
   // int lastNumber = 300;
    double start_time, run_time;

    start_time = omp_get_wtime();
    int found = 0;

    // Initialize
    char* isPrime = new char[lastNumber + 1];
    for (int i = 0; i <= lastNumber; i++)
        isPrime[i] = 1;

    // Find all non-primes
    for (int i = 2; i * i <= lastNumber; i++)
        if (isPrime[i])
            for (int j = i * i; j <= lastNumber; j += i)
                isPrime[j] = 0;

    // Count prime numbers
    for (int i = 2; i <= lastNumber; i++) {

        found += isPrime[i];
    }



    cout << found;
    run_time = omp_get_wtime() - start_time;
    cout << "\nSecond : " << run_time << "s";

    delete[] isPrime;
    return found;
}


int main()
{
    //int lastNumber = 300;
    double start_time, run_time;

    start_time = omp_get_wtime();
    int found = 0;

    // Initialize
    char* isPrime = new char[lastNumber + 1];
    
    #pragma omp parallel for
    for (int i = 0; i <= lastNumber; i++)
        isPrime[i] = 1;
    
    // Find all non-primes
    #pragma omp parallel for schedule(dynamic) collapse(2)
    for (int i = 2; i  <= (int)sqrt((double)lastNumber); i++)
        if (isPrime[i])
            for (int j = i * i; j <= lastNumber; j += i)
                isPrime[j] = 0;

    // Count prime numbers
   #pragma omp parallel for reduction(+:found)
    for (int i = 2; i <= lastNumber; i++) {

        found += isPrime[i];
    }

    cout << "There are total of " << found << " prime number found within " << lastNumber;
    run_time = omp_get_wtime() - start_time;
    cout << "\nSecond : " << run_time << "s";

    delete[] isPrime;
    return found;
}

int main3()
{
    int i, nthreads;
    double start_time, run_time;
    omp_set_num_threads(NUM_THREADS);

    start_time = omp_get_wtime();
    int found = 0;

    // Initialize
    char* isPrime = new char[lastNumber + 1];

    #pragma omp parallel
    {
        int i, id, nthrds;
        id = omp_get_thread_num(); // based on parallel
        nthrds = omp_get_num_threads(); // 4

        for (i = id; i <= lastNumber; i = i + nthrds)
            isPrime[i] = 1;
    }

   
    // Find all non-primes
    #pragma omp parallel
    {
        int i, id, nthrds;
        id = omp_get_thread_num(); //based on parallel
        nthrds = omp_get_num_threads(); //4

        for (i = id;  i <= (int)sqrt((double)lastNumber); i = i + nthrds)
        {
            if (i == 1 || i == 0)
                continue;
            if (isPrime[i])
                for (int j = i * i; j <= lastNumber; j += i)
                    isPrime[j] = 0;
        }
    }


    // Count prime numbers
    #pragma omp parallel
    {
        int i, id, nthrds;
        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();     
        
        for (i = id; i <= lastNumber; i = i + nthrds) {
           // if (i == 1 || i == 0)
           //     continue;
            #pragma omp atomic
            found += isPrime[i];
        }
    }


    cout << "There are total of " << found-2 << " prime number found wiithin " << lastNumber;
    run_time = omp_get_wtime() - start_time;
    cout << "\nSecond : " << run_time << "s";

    delete[] isPrime;
    return found;
}


int main4(int argc, char** argv) {
    /* Declare variables */
    int N = 300000000; /* The positive integer under which we are finding primes */
    int sqrtN = 0; /* The square root of N, which is stored in a variable to
                      avoid making excessive calls to sqrt(N) */
    int c = 0; /* Used to check the next number to be circled */
    int m = 0; /* Used to check the next number to be marked */
    char* list1; /* The list of numbers <= sqrtN -- if list1[x] equals 1, then x
                   is marked.  If list1[x] equals 0, then x is unmarked. */
    char* list2; /* The list of numbers > sqrtN -- if list2[x-L] is marked, then
                   x is marked.  If list2[x-L] equals 0, then x is unmarked. */
    char next_option = ' '; /* Used for parsing command line arguments */
    int S = 0; /* A near-as-possible even split of the count of numbers above
                  sqrtN */
    int R = 0; /* The remainder of the near-as-possible even split */
    int L = 0; /* The lowest number in the current process's split */
    int H = 0; /* The highest number in the current process's split */
    int r = 0; /* The rank of the current process */
    int p = 0; /* The total number of processes */
    double start_time, run_time;
    int found = 0, sum = 0; int error;
    /* Initialize the MPI Environment */
    MPI_Init(&argc, &argv);

    /* Determine the rank of the current process and the number of processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Comm_size(MPI_COMM_WORLD, &p);


    /* Calculate sqrtN */
    sqrtN = (int)sqrt(N);

    /* Calculate S, R, L, and H */

    // Split size for different processes 
    S = (N - (sqrtN + 1)) / p;

    // Remainder size for last processes
    R = (N - (sqrtN + 1)) % p;

    // Lowest boundary in each split
    L = sqrtN + r * S + 1;

    // Highest boundary in each split
    H = L + S - 1;

    // finally until the last process, responsible for the remainder
    if (r == p - 1) {
        H += R;
    }

    /* Allocate memory for lists */
    list1 = (char*)malloc((sqrtN + 1) * sizeof(char));
    list2 = (char*)malloc((H - L + 1) * sizeof(char));

    /* Exit if malloc failed */
    if (list1 == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for list1.\n");
        exit(-1);
    }
    if (list2 == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for list2.\n");
        exit(-1);
    }

    start_time = omp_get_wtime();

    //initialize array with 1 
    // assume all is one first 
    //list1 is	designated	for	the	numbers 2 through sqrtN
    for (c = 2; c <= sqrtN; c++) {
        list1[c] = 1;
    }

    //list2 contains a section	of	the	remaining	numbers. 
    // different processes will run the program once, each has different L and H
    for (c = L; c <= H; c++) {
        list2[c - L] = 1;
    }

    /* Run through each number in list1 */
    for (c = 2; c <= sqrtN; c++) {
        /* If the number is marked */
        if (list1[c]) {
            /* Run through each number bigger than c in list1 */
            for (m = c + c; m <= sqrtN; m += c) {
                /* If m is a multiple of c */
                    /* Mark m */
                list1[m] = 0;
            }
            //===========================================
            //=============================================
            int temp = L;
            /* Run through each number bigger than c in list2 */
            while (temp % c != 0) {
                temp++;
            }

            for (m = temp; m <= H; m += c)
            {
                // m = L = 4 // m - L = 0
                // m = 6 L = 4 // c = 2 // m- L = 2
                list2[m - L] = 0;
            }
        }
    }

    
    /* If Rank 0 is the current process */
    if (r == 0) {
        /* Run through each of the numbers in list1 */
        for (c = 2; c <= sqrtN; c++) {
            /* If the number is unmarked */
            if (list1[c] == 1) {
                /* The number is prime, print it */
               // printf("%d ", c);
                found++;
            }
        }

        /* Run through each of the numbers in list2 */
        for (c = L; c <= H; c++) {
            /* If the number is unmarked */
            if (list2[c - L] == 1) {
                /* The number is prime, print it */
               // printf("%d ", c);
                found++;
            }
        }

        /* Run through each of the other processes */
        for (r = 1; r <= p - 1; r++) {

            /* Calculate L and H for r */
            L = sqrtN + r * S + 1;
            H = L + S - 1;

            // assigtn to last processes
            if (r == p - 1) {
                H += R;
            }

            /* Receive list2 from the process */
            MPI_Recv(list2, H - L + 1, MPI_CHAR, r, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);

            /* Run through the list2 that was just received */
            for (c = L; c <= H; c++) {

                /* If the number is unmarked */
                if (list2[c - L] == 1) {
                    /* The number is prime, print it */
                   // printf("%d ", c);
                    found++;
                }
            }
        }
        
        
        MPI_Send(&found, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        printf("\n");

        /* If the process is not Rank 0 */
    }
    else {

        /* Send list2 to Rank 0 */
        MPI_Send(list2, H - L + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }

    run_time = omp_get_wtime() - start_time;

    MPI_Barrier(MPI_COMM_WORLD);
    if (r == 1) {

        MPI_Recv(&found,  1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("\n\nTotal elapsed time: %10.6f seconds\n", run_time);
        cout << "There are total of " << found << " prime number found within " << N << "\n\n";
    }
    /* Deallocate memory for list */
    if (r == 0) {
        free(list2);
        free(list1);
    }
    
    //       MPI_Barrier(MPI_COMM_WORLD);
    /* Finalize the MPI environment */
    error = MPI_Finalize();

    return 0;
}
