
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

// Makes program global
unsigned int seed = 0;

// Global State Variables
int	NowYear;		// 2024- 2029
int	NowMonth;		// 0 - 11

float	NowPrecip;		// inches of rain per month
float	NowTemp;		// temperature this month
float	NowHeight;		// grain height in inches
int	NowNumDeer;		// number of deer in the current population
int     NowNumBear;            // number of bear in the current population

// Basic time step of one month paramaters:
const float GRAIN_GROWS_PER_MONTH =	       12.0;
const float ONE_DEER_EATS_PER_MONTH =		1.0;
const int ONE_BEAR_EATS_DEER_PER_MONTH =          2;
const float ONE_BEAR_EATS_GRAIN_PER_MONTH =     1.0;

const float AVG_PRECIP_PER_MONTH =		7.0;	// average
const float AMP_PRECIP_PER_MONTH =		6.0;	// plus or minus
const float RANDOM_PRECIP =			2.0;	// plus or minus noise

const float AVG_TEMP =				60.0;	// average
const float AMP_TEMP =				20.0;	// plus or minus
const float RANDOM_TEMP =			10.0;	// plus or minus noise

const float MIDTEMP =				40.0;
const float MIDPRECIP =				10.0;

// Barrier solution
omp_lock_t	Lock;
volatile int	NumInThreadTeam;
volatile int	NumAtBarrier;
volatile int	NumGone;

// Random number generator using Ranf() function from proj01
float
Ranf( float low, float high )
{
        float r = (float) rand();               // 0 - RAND_MAX
        float t = r  /  (float) RAND_MAX;       // 0. - 1.

        return   low  +  t * ( high - low );
}

void
TimeOfDaySeed( )
{
	struct tm y2k = { 0 };
	y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
	y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;

	time_t  timer;
	time( &timer );
	double seconds = difftime( timer, mktime(&y2k) );
	unsigned int seed = (unsigned int)( 1000.*seconds );    // milliseconds
	srand( seed );
}

// Returns square root
float
SQR( float x )
{
        return x*x;
}

// Function prototypes
void	InitBarrier( int );
void	WaitBarrier( );

// Specify how many threads will be in the barrier:
//	(also init's the Lock)
void
InitBarrier( int n )
{
        NumInThreadTeam = n;
        NumAtBarrier = 0;
	omp_init_lock( &Lock );
}

// Makes calling thread wait until all other threads have caught up:
void
WaitBarrier( )
{
        omp_set_lock( &Lock );
        {
                NumAtBarrier++;
                if( NumAtBarrier == NumInThreadTeam )
                {
                        NumGone = 0;
                        NumAtBarrier = 0;
                        // let all other threads get back to what they were doing
			// before this one unlocks, knowing that they might immediately
			// call WaitBarrier( ) again:
                        while( NumGone != NumInThreadTeam-1 );
                        omp_unset_lock( &Lock );
                        return;
                }
        }
        omp_unset_lock( &Lock );

        while( NumAtBarrier != 0 );	// this waits for the nth thread to arrive

        #pragma omp atomic
        NumGone++;			// this flags how many threads have returned
}

// Function for temp and precipitation of particular month
void
TempPrecip()
{
        float ang = (  30.*(float)NowMonth + 15.  ) * ( M_PI / 180. );	// angle of earth around the sun

        float temp = AVG_TEMP - AMP_TEMP * cos( ang );
        NowTemp = temp + Ranf( -RANDOM_TEMP, RANDOM_TEMP );

        float precip = AVG_PRECIP_PER_MONTH + AMP_PRECIP_PER_MONTH * sin( ang );
        NowPrecip = precip + Ranf( -RANDOM_PRECIP, RANDOM_PRECIP );
        if( NowPrecip < 0. )
	        NowPrecip = 0.;
}

// Parallel Sections Directive
void
Deer()
{
        while( NowYear < 2030 )
        {
	        // compute a temporary next-value for this quantity
	        // based on the current state of the simulation:
	        int nextNumDeer = NowNumDeer - ( NowNumBear * ONE_BEAR_EATS_DEER_PER_MONTH );
                int carryingCapacity = (int)( NowHeight );
                
                // Calculating next number of deer
                if( nextNumDeer < carryingCapacity )
                        nextNumDeer++;
                else if( nextNumDeer > carryingCapacity )
                {
                        nextNumDeer--;
                }        

                if( nextNumDeer < 0 )
                {
                        nextNumDeer = 0;
                }

	        // DoneComputing barrier:
	        WaitBarrier( );
	        NowNumDeer = nextNumDeer;

	        // DoneAssigning barrier:
	        WaitBarrier( );

	        // DonePrinting barrier:
	        WaitBarrier( );
	        
        }
}

void
Grain()
{
        while( NowYear < 2030)
        {
               // compute a temporary next-value for this quantity
	        // based on the current state of the simulation:
                float tempFactor = exp(   -SQR(  ( NowTemp - MIDTEMP ) / 10.  )   );
                float precipFactor = exp(   -SQR(  ( NowPrecip - MIDPRECIP ) / 10.  )   );
                float nextHeight = NowHeight;
                
                // Calculating next height of grain
                nextHeight += tempFactor * precipFactor * GRAIN_GROWS_PER_MONTH;
                nextHeight -= (float)NowNumDeer * ONE_DEER_EATS_PER_MONTH;
                nextHeight -= (float)NowNumBear * ONE_BEAR_EATS_GRAIN_PER_MONTH;

                if ( nextHeight <= 0. )
                {
                        nextHeight = 0.;
                }

	        // DoneComputing barrier:
	        WaitBarrier( );
	        NowHeight = nextHeight;

	        // DoneAssigning barrier:
	        WaitBarrier( );

	        // DonePrinting barrier:
	        WaitBarrier( );
        }
}

void
Watcher()
{
        while( NowYear < 2030)
        {
                // compute a temporary next-value for this quantity
	        // based on the current state of the simulation:

                // DoneComputing barrier:
	        WaitBarrier( );

	        // DoneAssigning barrier:
	        WaitBarrier( );

                fprintf(stderr, "%2d-%2d , %6.f , %6.f , %6.f , %6d , %2d\n",
                        NowMonth + 1, NowYear, NowPrecip, NowTemp, NowHeight, NowNumDeer, NowNumBear);
                
                // Calculating month parameter
                if (NowMonth >= 11)
                {
                        NowMonth = 0;
                        NowYear++;
                }
                else
                {
                        NowMonth++;
                }

                // Calculating Environmental parameters
                TempPrecip( );
                TimeOfDaySeed( );

	        // DonePrinting barrier:
	        WaitBarrier( );
        }
}

void
Bear()
{
        while( NowYear < 2030)
        {
                // compute a temporary next-value for this quantity
	        // based on the current state of the simulation:
                int nextNumBear = NowNumBear;
                int carryingCapacity = NowNumDeer;

                // Calculating next number of bears 
                if( nextNumBear < carryingCapacity )
                {
                        nextNumBear++;
                }      
                else if( nextNumBear > carryingCapacity )
                {
                        nextNumBear--;
                }              

                if( nextNumBear < 0 )
                {       
                        nextNumBear = 0;
                }

	        // DoneComputing barrier:
	        WaitBarrier( );
                NowNumBear = nextNumBear;

	        // DoneAssigning barrier:
	        WaitBarrier( );

	        // DonePrinting barrier:
	        WaitBarrier( );
        }
}

// Main program
int
main( int argc, char *argv[ ] )
{
#ifdef _OPENMP
	#ifndef CSV
		fprintf( stderr, "OpenMP is supported -- version = %d\n", _OPENMP );
	#endif
#else
        fprintf( stderr, "No OpenMP support!\n" );
        return 1;
#endif

        // starting date and time:
        NowMonth =    0;
        NowYear  = 2024;

        // starting state (feel free to change this if you want):
        NowNumDeer = 2;
        NowNumBear = 1;
        NowHeight =  5.;
        
        // Set Environmental Parameters
        TempPrecip( );
        TimeOfDaySeed( );

        omp_set_num_threads( 4 );	// same as # of sections
        InitBarrier( 4 );
        #pragma omp parallel sections
        {
	        #pragma omp section
	        {
	        	Deer( );
	        }

	        #pragma omp section
	        {
	        	Grain( );
	        }

	        #pragma omp section
	        {
	        	Watcher( );
	        }

	        #pragma omp section
	        {
	        	Bear( );	// your own
	        }
        }       // implied barrier -- all functions must return in order
	        // to allow any of them to get past here
}