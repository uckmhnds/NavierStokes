#ifndef READMESH_H_INCLUDED
#define READMESH_H_INCLUDED

void readMesh(){

    std::ifstream gridFile;
    gridFile.open(filename);

    if (!gridFile.is_open()){
        std::cout << "Mesh file is not found..." << std::endl;
    }else{
        std::cout << "Reading mesh file..." << std::endl;
    }

    gridFile >> nodeNumber >> cellNumber;

    for (int i=0; i<nodeNumber; i++){
        gridFile >> nodes[i].x >> nodes[i].y;
    }

    for (int i=0; i<cellNumber; i++){
        gridFile >> cells[i].nodeIDs[0] >> cells[i].nodeIDs[1] >> cells[i].nodeIDs[2] >> cells[i].neighID[0] >> cells[i].neighID[1] >> cells[i].neighID[2];
    }

    float dx, dx1;
    float dy, dy1;
    int j2;

    for (int i=0; i<cellNumber; i++){

        for (int j=0; j<3; j++){

            j2                      = (j + 1) % 3;

            dx                      = nodes[cells[i].nodeIDs[j2] - 1].x - nodes[cells[i].nodeIDs[j] - 1].x;
            dy                      = nodes[cells[i].nodeIDs[j2] - 1].y - nodes[cells[i].nodeIDs[j] - 1].y;

            cells[i].surfaceN[j]    = sqrt(dx * dx + dy * dy);
            cells[i].surfaceC[j]    =  dy / cells[i].surfaceN[j];
            cells[i].surfaceS[j]    = -dx / cells[i].surfaceN[j];

        }

        cells[i].centerX            = (nodes[cells[i].nodeIDs[0] - 1].x + nodes[cells[i].nodeIDs[1] - 1].x + nodes[cells[i].nodeIDs[2] - 1].x) / 3.;
        cells[i].centerY            = (nodes[cells[i].nodeIDs[0] - 1].y + nodes[cells[i].nodeIDs[1] - 1].y + nodes[cells[i].nodeIDs[2] - 1].y) / 3.;

        dx1                         = nodes[cells[i].nodeIDs[1] - 1].x - nodes[cells[i].nodeIDs[0] - 1].x;
        dy1                         = nodes[cells[i].nodeIDs[1] - 1].y - nodes[cells[i].nodeIDs[0] - 1].y;

        cells[i].area               = 0.5 * (dx * dy1 - dx1 * dy);
        cells[i].dxmin              = 0.7 * std::min(std::min(cells[i].surfaceN[0], cells[i].surfaceN[1]), cells[i].surfaceN[2]);

        if (cells[i].area < 0.){std::cout << "*\n**\n***\n****\nZERO VOLUME\n****\n***\n**\n*\n";}

    }

}


#endif // READMESH_H_INCLUDED
