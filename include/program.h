#include <string>
using namespace std;

#ifndef PROGRAM_H
#define PROGRAM_H


class program
{
    public:
        program();
        virtual ~program();
        void run (char* xml_input);
        void copyfile( char* SRC, std::string  DEST);
    protected:
    private:
};

#endif // PROGRAM_H
