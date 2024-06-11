#include <iostream>
#include <math.h>
#include <cstdlib> // for rand() and srand()
#include <ctime>   // for time()
#include <iomanip> // for setprecision


#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))
#define Sqr(xx)     ((xx) * (xx))
#define Cube(xx)    ((xx) * (xx) * (xx))
#define AllocMem(a, n, t)  a = (t *) malloc ((n) * sizeof (t))


using namespace std;

const int n = 100;   // Number of particles
const double l = 10.0;   // Length of the box

//Random number generator
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        float temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        idum=IA1(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1;
        idum=IA1(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}



// Function to initialize particle positions within the box in a lattice
void initializeParticles(double p[][3], int size, double length) {
    double spacing = length / cbrt(size);
    int count = 0;

    for (double x = 0; x < length && count < size; x += spacing) {
        for (double y = 0; y < length && count < size; y += spacing) {
            for (double z = 0; z < length && count < size; z += spacing) {
                if (count < size) {
                    p[count][0] = x;
                    p[count][1] = y;
                    p[count][2] = z;
                    ++count;
                }
            }
        }
    }
}

// Function to compute energy
double computeEnergy(double p[][3], int size) {
    double epsilon = 1.0;
    double r;
    double t;
    int sigma = 1;
    double rc = (2.5 * sigma);
    double E = 0;
    double E_cutoff;
    double y = pow((sigma / rc), 6);
    double dx,dy,dz;
    E_cutoff = 4 * epsilon * ((y * y) - y);
    for (int i = 0; i < (size - 1); ++i) {
        for (int j = (i + 1); j < size; ++j) {
            dx=p[j][0] - p[i][0];
            dy=p[j][1] - p[i][1];
            dz=p[j][2] - p[i][2];
            dx-=l*NINT(dx/l);
            dy-=l*NINT(dy/l);
            dz-=l*NINT(dz/l);
            r = Sqr(dx) + Sqr(dy) + Sqr(dz);
            t = pow((sigma * sigma / r), 3);
            E += (4 * epsilon * ((t * t) - t)) - E_cutoff;
        }
    }
    return E;
}

// Function to displace the particle slightly
void updatePositions(double p[][3], int id, double length) {
    double delta;
    for (int j = 0; j < 3; ++j) {
        delta = (rand() / (double)RAND_MAX) - 0.5;  // Range from -0.5 to 0.5
        p[id][j] += delta;
        if (p[id][j] > length) {
            p[id][j]-= length;
        }
        if (p[id][j] < 0.0) {
            p[id][j]+= length;
        }
    }
}

// Function to display initial particle positions
void initialParticlesPosition(double p[][3], int size, double length) {
    cout << size << endl;
    cout << "0" << endl;
    for (int i = 0; i < size; ++i) {
        cout << "0" << "\t";
        for (int j = 0; j < 3; ++j) {
            cout << fixed << setprecision(4) << p[i][j] << "\t\t";
        }
        cout << "\n";
    }
}

// Function to write the position of the particles after updating them
void write(double p[][3], int k, int size) {
    cout << size << endl;
    cout << (k + 1) << endl;
    for (int i = 0; i < size; ++i) {
        cout << "0" << "\t";
        for (int j = 0; j < 3; ++j) {
            cout << fixed << setprecision(4) << p[i][j] << "\t\t";
        }
        cout << "\n";
    }
}

int main() {
    srand(time(0));
    int id;
    double p[n][3]; // Array to store particle positions
    double xo,yo,zo;
    FILE *ip,*op;
    long idum;
    idum=-12412424+rand()%100;
    initializeParticles(p, n, l);
    initialParticlesPosition(p, n, l);

    op=fopen("pos.xyz","w");
    fclose(op);


    int nstep = 1000; // Number of steps
    for (int i = 0; i < nstep; ++i) {
        for (int j = 0; j < n; ++j) {
            //id = (rand() / (double)RAND_MAX) * n;
            id = (int)(ran2(&idum)*(double)n);
            double e1 = computeEnergy(p, n);
            xo=p[id][0];
            yo=p[id][1];
            zo=p[id][2];
            updatePositions(p, id, l);
            double e2 = computeEnergy(p, n);
            double d = e2 - e1;
            //if (d < 0 || (exp(-d) > (rand() / (double)RAND_MAX))) {
            //}
            //else
            if(d>0)
            if( ran2(&idum)>exp(-d))
            {
                p[id][0]=xo;
                p[id][1]=yo;
                p[id][2]=zo;
            }
        }
        if (i % 10 == 0) {
            //write(p, i, n);
        printf("%d\n",i);
        op=fopen("pos_eps10e-1.xyz","a");
        fprintf(op,"%d\n%d\n",n,i);
            for(int k=0;k<n;k++)
                fprintf(op,"1 %lf %lf %lf\n",p[k][0],p[k][1],p[k][2]);
        fclose(op);

        }
    }
    return 0;
}
