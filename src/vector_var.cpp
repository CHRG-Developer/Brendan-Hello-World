#include "vector_var.h"
#include <math.h>       /* pow ,cosh */

vector_var::vector_var()
{
    //ctor
}

//vector_var::~vector_var()
//{
//    //dtor
//}

double  vector_var::Dot_Product(vector_var b){

    return x * b.x + y *b.y + z * b.z;

}
double vector_var::Magnitude(){

    return sqrt(pow(x,2) + pow(y,2) + pow(z,2) );
}
double  vector_var::Angle_Between_Vectors(vector_var b){
    double num, denom;
    num = Dot_Product(b);
    denom = Magnitude();
    denom = denom* b.Magnitude();
    num = num/denom;
    num = acos (num);

    return acos (Dot_Product(b)/(Magnitude()* b.Magnitude())) ;

}


void vector_var::Get_Gradient(double y1, double y2, vector_var x1, vector_var x2 ){

        x = y2 -y1 ;
        y = y2 -y1 ;

        //check if x1 and x2 are same co-ordinate -> zero gradient
        if( fabs(x2.x - x1.x ) < pow(10,-5) ){
            x = 0;
        }else{
            x =  x / (x2.x - x1.x );
        }
        if( fabs(x2.y -x1.y) < pow(10,-5) ){
            y = 0;
        }else{
            y = y / (x2.y -x1.y);
        }
        z = 0; //update for 3d


}
