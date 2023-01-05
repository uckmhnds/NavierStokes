#ifndef OBJECTS_H_INCLUDED
#define OBJECTS_H_INCLUDED

//
class Cell{

    public:
        int nodeIDs[3];
        int neighID[3];
        float surfaceN[3];
        float surfaceC[3];
        float surfaceS[3];
        float centerX;
        float centerY;
        float area;
        float dxmin;
        Cell () {};

};
//
class Node{

    public:
        float x, y;
        Node(){};
};
//
class Q{

    public:
        float rho;
        float rhou;
        float rhov;
        float e;
        float p;
        float T;
};
//
class Flux{

    public:
        float f1;
        float f2;
        float f3;
        float f4;
};

Node * nodes    = new Node[nodeNumber];
Cell * cells    = new Cell[cellNumber];
Q    * q        = new Q[cellNumber];
Q    * q_       = new Q[cellNumber];
Q    * dq       = new Q[cellNumber];
Q    * qnode    = new Q[nodeNumber];
Q    * qnode_   = new Q[nodeNumber];
Flux * tflux    = new Flux[cellNumber];
Flux * tflux_   = new Flux[cellNumber];

Q ql_, qr_;
Flux flux_, flux__;

#endif // OBJECTS_H_INCLUDED
