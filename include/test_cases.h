#ifndef TEST_CASES_H
#define TEST_CASES_H


class test_cases
{
    public:
        test_cases();
        virtual ~test_cases();
        void west_to_east_couette_flow();
        void east_to_west_couette_flow();
        void north_to_south_couette_flow();
        void south_to_north_couette_flow();
        void vector_var_tests();
        void lid_driven_cavity_N();
    protected:
    private:
        double X,Y,dx,dy,dt;
        double reynolds,kine_viscosity;
        double U;
        double simulation_length;

};

#endif // TEST_CASES_H
