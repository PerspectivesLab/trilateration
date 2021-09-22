// trilateration.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
#include <math.h> 
#include <cmath>


class vector3 {

public:


    vector3() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }


    vector3( double xIn, double yIn, double zIn ) { 

        x = xIn;
        y = yIn;
        z = zIn;
    }

    double x;
    double y;
    double z;


    double r;

};




double dot(vector3 a, vector3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vector3   vector_subtract(vector3 a, vector3 b)
{
    return   vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}


vector3 vector_add(vector3 a, vector3 b)
{
    return   vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}


vector3 vector_divide(vector3 a, double b)
{
    return   vector3(a.x / b, a.y / b, a.z / b);
}


vector3 vector_multiply(vector3 a, double b)
{
    return   vector3(a.x * b, a.y * b, a.z * b);
}


vector3 vector_cross(vector3 a, vector3 b)
{
    return   vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z , a.x * b.y - a.y * b.x);
}


double sqr(double a)
{
    return   a* a;
}

double norm(vector3 a)
{
    return sqrt(sqr(a.x) + sqr(a.y) + sqr(a.z));
}

 
void trilaterate(vector3 p1, vector3 p2, vector3 p3, vector3 & res1, vector3& res2)
{
 
    vector3 ex, ey, ez, a;// , j, , a, x, y, z, b, p4;

    ex = vector_divide(vector_subtract(p2, p1), norm(vector_subtract(p2, p1)));

    double i = dot(ex, vector_subtract(p3, p1));
    a = vector_subtract(vector_subtract(p3, p1), vector_multiply(ex, i));
    ey = vector_divide(a, norm(a));
    ez = vector_cross(ex, ey);
    double d = norm(vector_subtract(p2, p1));
    double j = dot(ey, vector_subtract(p3, p1));

    double x = (sqr(p1.r) - sqr(p2.r) + sqr(d)) / (2 * d);
    double y = (sqr(p1.r) - sqr(p3.r) + sqr(i) + sqr(j)) / (2 * j) - (i / j) * x;

    double b = sqr(p1.r) - sqr(x) - sqr(y);

    // floating point math flaw in IEEE 754 standard
    // see https://github.com/gheja/trilateration.js/issues/2
    if (abs(b) < 0.0000000001)
    {
        b = 0;
    }

    double z = sqrt(b);

    // no solution found
    if (isnan(z))
    {

        throw std::invalid_argument("no solution found");
        
    }

    vector3 res  = vector_add(p1, vector_add(vector_multiply(ex, x), vector_multiply(ey, y)));
    vector3 p4a = vector_add(res, vector_multiply(ez, z));
    vector3 p4b = vector_subtract(res, vector_multiply(ez, z));

    if (z == 0 ) 
    {
        res1 = res;
    }
    else
    {
        // 2 solution par rapport a la face du triangle
        res1 = p4a;
        res2 = p4b;
        
    }
}



int main()
{
    std::cout << "Hello World!\n";

    // position en Z
    double sensorZ = 10;
    
    
    vector3 p1(1.0, 1.0, sensorZ);
    p1.r = 4.042;// distance metrique du signal par rapport au point

    vector3 p2(2.0, 2.0, sensorZ);
    p2.r = 4.046;

    vector3 p3(3.0, 1.0, sensorZ);
    p3.r = 4.01;

    

    vector3  res1;
    vector3  res2;
    // !! les points doivent se suivre dans le sens anti horaire !!
    trilaterate(p1, p3, p2, res1, res2);

      
    std::cout << "estimation de la position 3D: " << res1.x << " " << res1.y << " " << res1.z;


    getchar();


}

 