#ifndef VECTOR_VAR_H
#define VECTOR_VAR_H


class vector_var
{
    public:
        vector_var();
        virtual ~vector_var(){};
        double x;
        double y;
        double z;
        double Dot_Product(vector_var b);
        double Magnitude();
        double Angle_Between_Vectors( vector_var b);
        void Get_Gradient(double y1, double y2, vector_var x1, vector_var x2 );
    protected:
    private:
};

#endif // VECTOR_VAR_H
