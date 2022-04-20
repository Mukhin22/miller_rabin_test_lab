#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

/* How many numbers is needed to find with Miller Rabin test*/
#define PRIME_NUMS_TO_FIND  5
/* My variant input parameters*/
#define NUM_TO_START_FROM   ((int64_t) 595000000)
#define ACCURACY            ((double)  0.0001)

/* @brief Macro is needed to build the coef of RAND_MAX
 * Linux GCC compiler has higher value and gives
 * correct results but windows minGW has RAND_MAX
 * about 32k - so result was incorrect
 * @retval - 1 in case RANMAX is high enough
 * and double coef in case if it's needed (minGW)
 */
#define SYS_RAND_MAX_COEF(_num)  ((RAND_MAX < (_num)) ? \
  ((double)(((double)_num) / ((double)RAND_MAX))) : 1 )

// type to carry the launch parameters
typedef struct params {
	double accuracy;
	int64_t num_to_start_from;
} algo_params_t;

// start randomizer dependent on system time
static void rand_init(void)
{
	srand(time(0));
}

static int64_t round_to_higher(const double num)
{
	int64_t inum = (int64_t) num;
    if (num == (double)inum) {
        return inum;
    }
    return inum + 1;
}

static int64_t calc_tests_num(const double epsilon)
{
	// Math: log of base b on a = log a / log b
	double res = abs(log(epsilon) / log(4));
	return round_to_higher(res);
}

// exponentional modulo function. Formula: (base ^ exp) % mod
static int64_t powmod(int64_t base, int64_t exp, int64_t mod)
{
	int64_t res = 1;

	while (exp != 0) {
		if ((exp & 1) != 0) {
			res = (1ll * res * base) % mod;
		}

		base = (1ll * base * base) % mod;
		exp >>= 1;
	}
	
	return res;
}
static inline int64_t find_d(const int64_t n)
{
	// Find r such that n = 2^d * r + 1 for some r >= 1
    int64_t d = n - 1;
    while (d % 2 == 0) {
        d /= 2;
    }
    return d;
}

static int64_t n0 = 0; // used for tests of display console to form tables
 					   // printf("%ld\n", a); this one is for print when used

// main miller rabin test function checking for conditions (1) and (2) for lab
// use true for display flag if printing of numbers checked is needed to console
static bool miller_rabin_test(int64_t d, int64_t n, bool display, FILE *dump_file)
{
    // Pick a random number in [2..n-1]
    int64_t cut_case = 0;
    int64_t a = 2 + ((rand() % (n - 1)) * SYS_RAND_MAX_COEF(n - 1));

 	if (display && (n0 != n)) {
 		// n0 = n;
 		printf("Testing witness %ld for number %ld\n", a, n);
 		fprintf(dump_file, "Testing witness %ld for number %ld\n", a, n);
 	}
    // Compute a^d % n
    int64_t x = powmod(a, d, n);
 
    if ((x == 1) || (x == (n - 1))) {
    	if (display) {
    		printf("Number %ld is a witness of simplicity for number %ld\n", a, n);
    		fprintf(dump_file, "Number %ld is a witness of simplicity for number %ld\n", a, n);
 			// printf("%ld\n", a);
 		}
       return true;
    }
 
    // Keep squaring x while one of the following doesn't
    // happen
    // (a)   d does not reach n-1
    // (b)   (x^2) % n is not 1
    // (c)   (x^2) % n is not n-1
    while (d != (n - 1))
    {
        x = (x * x) % n; // fast powmod condition
        d *= 2;
 		
 		// cycle will always get next 1, so it's false condition
        if (x == 1) {
        	if (display) {
        		printf("Number %ld is a witness of compositness for number %ld\n", a, n);
        	    fprintf(dump_file, "Number %ld is a witness of compositness for number %ld\n", a, n);
        	}
        	return false;
        }
        // codition (2) happen 
        if (x == (n - 1)) {
        	if (display) {
				printf("Number %ld is a witness of simplicity for number %ld\n", a, n);
				fprintf(dump_file, "Number %ld is a witness of simplicity for number %ld\n", a, n);  			  			
			}
        	return true;
    	}
    }
    // Return composite
    if (display) {
    	printf("Number %ld is a witness of compositness for number %ld\n", a, n);
    	fprintf(dump_file, "Number %ld is a witness of compositness for number %ld\n", a, n);
    }
    return false;
}

// test the number according to lab conditions:
// number of tests, display info, etc.
static void is_prime_mil_rab_lab(int64_t tests_num, int64_t num_to_test,
	                             bool display, FILE *results_file,
	                             FILE *dump_file)
{
	int64_t d = 0;
	uint8_t numbers_found = 0;
	bool is_prime = false;

	// don't use odd numbers
	if ((num_to_test % 2) == 0) {
		num_to_test++;
	}

	// cycle to find needed count of prime numbers
	while (numbers_found < PRIME_NUMS_TO_FIND) {
		d = find_d(num_to_test);
		// cycle to test number simplicity
		for (int64_t i = 0; i < tests_num; ++i) {
				is_prime = miller_rabin_test(d, num_to_test, display, dump_file);
				if (!is_prime) {
					num_to_test += 2; // adding 2 to scip odd numbers
					d = find_d(num_to_test);
					i = -1;
					continue;
				}
		}
		// found one number let's try to find next
		numbers_found++;
		if (display) {
			printf("%ld value is prime\n", num_to_test);
			fprintf(dump_file, "%ld value is prime\n", num_to_test);
		}
		fprintf(results_file, "%ld\n", num_to_test); // write to file 
		num_to_test += 2;
	}
}

static inline void print_help(void)
{
	printf("To run the program correct use:\n"
		"./mil_rab.exe <int64_t number_to_start_from> <double accuracy>\n"
		"example: ./mil_rab.exe 50000000 0.00001\n");
}
static inline void print_error(void)
{
	printf("You've entered incorrect params, i will use defaults:"
	"number to start from  = %ld\n, accuracy = %lf\n", NUM_TO_START_FROM, ACCURACY);
}

// parser for start parameters
static void parse_params(int argc, char *argv[], algo_params_t* start_params)
{
	// check incorrect params number
	if (argc != 3) {
		print_error();
		start_params->num_to_start_from = NUM_TO_START_FROM;
		start_params->accuracy = ACCURACY;
		print_help();
		return;
	}

	int64_t start_num = 0;
	double accuracy = 0.0;
	accuracy = atof(argv[2]);
	start_num = atol(argv[1]);

	if (!accuracy || !start_num) {
		print_error();
		start_params->num_to_start_from = NUM_TO_START_FROM;
		start_params->accuracy = ACCURACY;
		print_help();
		return;		
	}

	start_params->num_to_start_from = start_num;
	start_params->accuracy = accuracy;
	printf("Parsed params are:\n"
		"accuracy = %lf, num_to_start_from = %ld\n", accuracy, start_num);
	printf("Parameters correct, let's find numbers\n");
}

int main(int argc, char *argv[])
{
	algo_params_t start_params = {};
	parse_params(argc, argv, &start_params);

	// File where final numbers stored
	FILE *results_file = fopen("results.txt", "w"); 
	if (!results_file) {
		printf("Error! Could not open results file\n"); 
		return -1;
	}

	//File where dump of logs is stored
	FILE *dump_file = fopen("mil_rab_cons_dump.txt", "w"); 
	if (!dump_file) {
		printf("Error! Could not open console dump file\n"); 
		return -1;
	}
	rand_init();

	int64_t tests_num = calc_tests_num(start_params.accuracy);
	printf("Number of tests to reach given accuracy is %ld \n\n", tests_num);
	fprintf(dump_file, "Number of tests to reach given accuracy is %ld \n\n", tests_num);

	is_prime_mil_rab_lab(tests_num, start_params.num_to_start_from, true,
	 					 results_file, dump_file);

	printf("Press any key to exit ...\n");
	getchar();

	fclose(results_file);
	fclose(dump_file);
	return 0;
}