#include "cgal_and_typedefs.h"
#include "logger.h"

int main() {
    logger.setLogging(true);
    logger << "Hello World!\n";
    logger << "A number: " << 2 << '\n';
}

