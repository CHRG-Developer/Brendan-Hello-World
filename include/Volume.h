#ifndef VOLUME_H
#define VOLUME_H


class Volume
{
    public:
        Volume();
        virtual ~Volume();

        double x,y,z;  // centroid of volume
        double get_x();
        double get_y();
        double get_z();


    protected:
    private:
};

#endif // VOLUME_H
