#ifndef INITIALIZE_H_INCLUDED
#define INITIALIZE_H_INCLUDED

void initialise(){

    std::cout << "Initializing..." << std::endl;

    for (int i=0; i<cellNumber; i++){

        q[i].rho    = rho_;
        q[i].rhou   = rho_ * u_;
        q[i].rhov   = rho_ * v_;
        q[i].e      = e_;
        q[i].p      = p_;
        q[i].T      = T_;

    }
}

#endif // INITIALIZE_H_INCLUDED
