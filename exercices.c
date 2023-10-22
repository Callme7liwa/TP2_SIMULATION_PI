/* 
    -- simPi : 
        Fonction pour estimer pi avec la méthode de Monte Carlo
  
    -- Entrees : nbPoints (type int)
       * @param nbPoints : le nombre de points utilisés pour l'estimation
    -- Sorties : estimation de pi de type double
         Renvoie l'estimation de pi  basée sur la méthode de Monte Carlo 
*/

/*   -- meanPi : 
        Fonction pour calculer la moyenne des estimations de Pi avec nombre definie d'experiences 
     
     -- Entrees : nbPoints (type int), nbExpriments (type int)
        * @param (nbPoints) Nombre de points de (type int) utilisé pour chaque estimation,
        * @param (nbExpriments) le nombre d'expériences (type int)  à effectuer.
     -- Sorties : 
        ->  Elle retourne la moyenne des estimations (type double) de Pi obtenues à partir des expériences.
*/

/*   -- meanPiWithResults : 
        Fonction pour calculer la moyenne des estimations de Pi avec nombre definie d'experiences toute en stockant 
        les valeurs de chaque simulation de PI dans un tableau resultat.
     
     -- Entrees : 
        * @param nbPoints         Le nombre de points générés dans chaque expérience.
        * @param nbExperiments    Le nombre total d'expériences effectuées.
     -- Sorties : 
        <- Elle retourne la moyenne des estimations (type double) de Pi obtenues à partir des expériences.
*/

/*   -- calculerMeanEtConfidenceRadius : 
        Calcule la moyenne de Pi et le rayon de confiance correspondant en utilisant les résultats d'un nombre donné d'expériences.

     -- Entrees : 
        * @param nbPoints         Le nombre de points générés dans chaque expérience.
        * @param nbExperiments    Le nombre total d'expériences effectuées.
        * @param ptrMeanPI        Un pointeur vers une variable où stocker la moyenne de Pi calculée.
        * @param ptrConfidenceRadius Un pointeur vers une variable où stocker le rayon de confiance calculé.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 

// Fonction pour calculer l'erreur absolue entre une estimation et la valeur réelle de pi
double absoluteError(double estimation)
{
    return fabs(estimation - M_PI);
}

// Fonction pour calculer l'erreur relative entre une estimation et la valeur réelle de pi
double relativeError(double estimation)
{
    return fabs(estimation - M_PI) / M_PI;
}

// Effectue une simulation pour estimer la valeur de pi en utilisant le nombre de points donné
double simPi(int nbPoints)
{
    int i, count = 0;
    double x, y;

    for (i = 0; i < nbPoints; i++)
    {
        x = genrand_real2();  // Génération d'un nombre aléatoire entre 0 et 1
        y = genrand_real2();

        if (x * x + y * y <= 1.0)
            count++;
    }

    return 4.0 * count / nbPoints;
}

// Effectue plusieurs simulations pour estimer la valeur de pi et retourne la moyenne des estimations
double meanPi(int nbPoints, int nbExperiments)
{
    double meanPI;
    double sommeCalcul = 0;

    for(int i = 0; i < nbExperiments; i++){
        sommeCalcul += simPi(nbPoints);
    }

    meanPI = sommeCalcul / (double)nbExperiments;

    return meanPI;
}

// Fonction pour calculer la moyenne des estimations de Pi avec nombre definie d'experiences toute en stockant 
// les valeurs de chaque simulation de PI dans un tableau resultat.
double meanPiWithResults(int nbPoints, int nbExperiments ,double resultats[])
{
    double meanPI;
    double sommeCalcul = 0;

    for(int i = 0; i < nbExperiments; i++){
        resultats[i] = simPi(nbPoints);
        sommeCalcul += resultats[i];
    }

    meanPI = sommeCalcul / (double)nbExperiments;

    return meanPI;
}

// Calcule la moyenne de Pi et le rayon de confiance correspondant en utilisant les résultats d'un nombre donné d'expériences.
void calculerMeanEtConfidenceRadius(int nbExperiments, int nbPoints, double *ptrMeanPI, double *ptrConfidenceRadius){
    
    double data[] = {
        12.706, 4.303, 3.182, 2.776, 2.571,
        2.447, 2.365, 2.308, 2.262, 2.228,
        2.201, 2.179, 2.160, 2.145, 2.131, 
        2.120, 2.110, 2.101, 2.093, 2.086,
        2.080, 2.074, 2.069, 2.064, 2.060, 
        2.056, 2.052, 2.048, 2.045, 2.042, 
        2.040, 2.038, 2.036, 2.034, 2.032,
        2.030, 2.028, 2.026, 2.024, 2.021,
        2.020, 2.018, 2.016, 2.014, 2.012,
        2.010, 2.008, 2.006, 2.004, 2.002,
        2.000, 2.000, 2.000, 2.000, 2.000, 
        2.000, 2.000, 2.000, 2.000, 2.000, 
        2.000, 2.000, 2.000, 2.000, 2.000, 
        2.000, 2.000, 2.000, 2.000, 2.000, 
        2.000, 2.000, 2.000, 2.000, 2.000, 
        2.000, 2.000, 2.000, 2.000, 2.000, 
        2.000, 2.000, 2.000, 2.000, 2.000, 
        2.000, 2.000, 2.000, 2.000, 2.000, 
        2.000, 2.000, 2.000, 2.000, 2.000, 
        2.000, 2.000, 2.000, 2.000, 2.000, 
        2.000, 2.000, 2.000, 2.000, 2.000, 
        2.000, 2.000, 2.000, 2.000, 2.000, 
        1.998, 1.996, 1.994, 1.992, 1.990,
        1.988, 1.986, 1.984, 1.982, 1.980,
        1.960
    };

    double *resultats = malloc(nbExperiments*sizeof(double));
    *ptrMeanPI = meanPiWithResults(nbPoints, nbExperiments, resultats);
    double estimateurSansBiais = 0;

    for(int i = 0; i < nbExperiments; i++){
        estimateurSansBiais += (resultats[i] - *ptrMeanPI) * (resultats[i] - *ptrMeanPI);
    }
    estimateurSansBiais = estimateurSansBiais / (nbExperiments-1);

    int dataSize = nbExperiments - 1;
    if(dataSize > 120){
        dataSize = 120;
    }

    *ptrConfidenceRadius = data[dataSize] * sqrt(estimateurSansBiais / nbExperiments);

    free(resultats);
}


// Cet exercice vise à estimer la valeur de Pi en utilisant la méthode de Monte Carlo avec différents niveaux de précision.
// Il utilise un générateur de nombres aléatoires et ajuste le nombre de points jusqu'à ce que la précision souhaitée soit atteinte.
void startExo1()
{
    int nbPoints[] = {1000, 1000000, 1000000000};
    int i, j;
    double estimationPi;
    
    for (i = 0; i < 3; i++)
    {
        int count = 0;
        double precision = 1.0;
        
        // Boucle pour répéter le calcul jusqu'à atteindre la précision souhaitée
        while (precision >= pow(10, -i-2))
        {
            // Appel de la fonction simPi pour estimer Pi avec le nombre de points actuel
            estimationPi = simPi(nbPoints[i]);
            
            // Calcul de la précision actuelle en comparant l'estimation de Pi avec la valeur réelle
            precision = fabs(estimationPi - M_PI);
            
            // Incrémentation du nombre de points
            nbPoints[i] += 1000;
            
            // Incrémentation du compteur de tentatives
            count++;
        }
        
        printf("Pour une précision de 10^-%d, il a fallu %d points pour estimer Pi à %lf\n", i+2, nbPoints[i], estimationPi);
    }
}

// L'objectif de la fonction startExo2 est de réaliser une série d'estimations de Pi en utilisant différents nombres de points 
// et d'expériences. Les résultats sont ensuite affichés, y compris la moyenne de Pi, l'erreur absolue et 
// l'erreur relative pour chaque combinaison de points et d'expériences.
void startExo2(){

    int nbPoints[] = {1000, 1000000, 1000000000};
    int nbExperiments[] = {10, 20, 30, 40};
    int i, j;

    for (i = 0; i < sizeof(nbPoints) / sizeof(nbPoints[0]); i++) {
        int currentNbPoints = nbPoints[i];
        printf("Estimations de Pi avec %d points :\n", currentNbPoints);
        for (j = 0; j < sizeof(nbExperiments) / sizeof(nbExperiments[0]); j++) {
            int currentNbExperiments = nbExperiments[j];
            double mean = meanPi(currentNbPoints, currentNbExperiments);
            double absError = absoluteError(mean);
            double relError = relativeError(mean);
            printf("Nombre d'expériences : %d\t Moyenne de Pi : %f\t Erreur absolue : %f\t Erreur relative : %f\n",
                    currentNbExperiments, mean, absError, relError);
        }
        printf("\n");
    }
}

// L'objectif de la fonction startExo3 est de calculer la moyenne de Pi et l'intervalle de confiance correspondant 
//pour un nombre donné d'expériences et de points. Les valeurs calculées sont ensuite affichées à l'écran. 
//Le code utilise la fonction calculerMeanEtConfidenceRadius pour effectuer le calcul de la moyenne et du rayon de confiance, puis les résultats sont affichés en utilisant printf.
void startExo3(){
    int nbExperiments = 40;
    int nbPoints = 1000;
    double meanPI,confidenceRadius;
    calculerMeanEtConfidenceRadius(nbExperiments, nbPoints, &meanPI, &confidenceRadius);
    printf("L'intervalle de confiance a 95 pourcents pour n = %d et nbrPoints = %d est de [%lf, %lf]\n", nbExperiments, nbPoints, meanPI - confidenceRadius, meanPI + confidenceRadius);
}


int main()
{
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    init_by_array(init, length);
    
    startExo1();
    // startExo2();
    // startExo3();
    return 0;
}