#include "Geometry.hpp"
#include "math.h"
#include "tol.hpp"
#include <iostream>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace FractureNetwork {

//calcola tutte le traccie e se sono passatnti o non passatni e le salva DOVE E COME??????
void CalculateTraces(vector<Fracture>& fractures, vector<Trace>& traces){
    unsigned int counter = 0; //tiene il conto tracce trovate (per ID)
    unsigned int dim = fractures.size(); //numero di fratture

    list<Trace> temp_traces = {}; // lista temporanea per salvare tracce
    vector<list<unsigned int>>  temp_passing(dim); // posizione i = lista id tracce passanti nel poligono i
    vector<list<unsigned int>>  temp_not_passing(dim); // posizione i = lista id tracce non passanti nel poligono i

    for(unsigned int i = 0; i<dim; i++){ //ciclo su tutte le matrici (poligoni)
        MatrixXd fracture1 = fractures[i].vertices;
        for(unsigned int j = i+1; j<dim; j++ ){ // ciclo su tutti i poligoni >i (altre coppie gia cocntrollate)
            MatrixXd fracture2 = fractures[j].vertices;
            if(near1(fracture1, fracture2)){ //se molto lontani non controllo nemmeno intersezione
                //inizializzo variabili utili
                Vector3d E1 = {};
                Vector3d E2 = {};
                bool tips1 = true;
                bool tips2 = true;
                bool find = fracturesIntersection(fracture1, fracture2, E1, E2,tips1, tips2); //calcolo intersezioni
                if(find){
                    //Salvo variabili in struttura adeguata
                    Trace traccia = {};

                    Matrix<double,3,2> vertices = {};
                    vertices.col(0) = E1;
                    vertices.col(1) = E2;
                    traccia.vertices = vertices;

                    traccia.ID = counter; //ID della traccia = numero in cui è stata trovata
                    traccia.first_generator = fractures[i].ID;
                    traccia.second_generator = fractures[j].ID;
                    double length = sqrt( pow((vertices(0,0)-vertices(0,1)),2) +  pow((vertices(1,0)-vertices(1,1)),2) +  pow((vertices(2,0)-vertices(2,1)),2)); //ATTENZIONE berrone no
                    traccia.length = length;

                    temp_traces.push_back(traccia);

                    //ATTENZIONE puschback diretto su vettore ridimensionando vettore???
                    if(!tips1){
                        temp_passing[i].push_back(counter);
                    }
                    else{
                        temp_not_passing[i].push_back(counter);
                    }
                    if(!tips2){
                        temp_passing[j].push_back(counter);
                    }
                    else{
                        temp_not_passing[j].push_back(counter);
                    }

                    counter ++;  
                }
            }
        }
    }

    //ATTENZIONE Vicini resize più sicuro ma si rompe con resize
    //ATTENZIONE serve sapere cosa fa move?
    traces.reserve(size(temp_traces));
    move(begin(temp_traces), end(temp_traces), back_inserter(traces));

    for(unsigned int i = 0; i<dim; i++){
        fractures[i].passing.reserve(temp_passing[i].size());
        move(begin(temp_passing[i]), end(temp_passing[i]), back_inserter(fractures[i].passing));
        fractures[i].not_passing.reserve(temp_not_passing[i].size());
        move(begin(temp_not_passing[i]), end(temp_not_passing[i]), back_inserter(fractures[i].not_passing));

    }
}

//***************************************************************************************************************
Vector3d normalP(const MatrixXd& fracture){
    Vector3d v1 = fracture.col(1)-fracture.col(0);
    Vector3d v2 = fracture.col(2)-fracture.col(1);
    return v1.cross(v2).normalized();
}

//***************************************************************************************************************
bool lineFractIntersect(const MatrixXd& fracture, const Vector3d& pointOnLine,const Vector3d& direction,vector<Vector3d>& intersections){

    for(unsigned int i=0;i<fracture.cols();i++){
        //intersezione lato-retta
        Vector3d ithFractureVertex = fracture.col(i);
        Vector3d ithFractureEdge =fracture.col((i+1)%fracture.cols())-fracture.col(i);

        if(abs(ithFractureEdge.cross(direction).norm())<tol2){
            if((ithFractureVertex-pointOnLine).cross(direction).isZero(tol2)){
                intersections.push_back(ithFractureVertex);
                intersections.push_back(fracture.col((i+1)%fracture.cols()));

            }
            continue;
        }
        double t = ((pointOnLine.cross(direction)).dot((ithFractureEdge.cross(direction)))-(ithFractureVertex.cross(direction)).dot((ithFractureEdge.cross(direction))))/(ithFractureEdge.cross(direction).norm()*ithFractureEdge.cross(direction).norm());
        if (t >= -tol2 && t <= 1.0+tol2){
            intersections.push_back(ithFractureVertex + t * ithFractureEdge);
        }
        if(intersections.size()==2){
            return true;
        }
    }
    return false;
}

//***************************************************************************************************************
bool fracturesIntersection(const MatrixXd& fracture1,const MatrixXd& fracture2, Vector3d& E1, Vector3d& E2, bool& tips1, bool& tips2)
{
    //Considero la retta di intersezione tra i due piani che contengono i poligoni
    //Calcolo la direction della retta come prodotto vettoriale tra le normali ai piani
    Vector3d n1 = normalP(fracture1);
    Vector3d n2 = normalP(fracture2);
    Vector3d direction = n1.cross(n2).normalized();
    //Se il prodotto vettoriale è nullo, i piani sono paralleli
    if(direction.isZero(tol1)){
        return false;
    }
    //Interseco il piano del poligono 1 con la retta che contiene il lato 1 del poligono 2 per trovare un punto della retta
    Vector3d edge = fracture2.col(1)-fracture2.col(0);
    Vector3d vertex = fracture2.col(0);
    //check che quel lato non sia parallelo al piano, e se lo è, prendo il lato successivo
  
    if(abs(edge.dot(n1))<tol1)
    {
        edge = fracture2.col(2)-fracture2.col(1);
        vertex = fracture2.col(1);
    }
    //Soluzione dell'intersezione piano-retta
    double t = -(n1.dot(vertex-fracture1.col(0)))/(n1.dot(edge));
    Vector3d pointOnLine = vertex+t*edge;
    //Ora controlliamo se il poligono interseca la retta

    vector<Vector3d> intersections;
    intersections.reserve(4);
    vector<Vector3d> inter1;
    if(lineFractIntersect(fracture1,pointOnLine,direction,inter1)){
        intersections.push_back(inter1[0]);
        intersections.push_back(inter1[1]);
    }
    else{
        return false;
    }
    vector<Vector3d> inter2;
    if(lineFractIntersect(fracture2,pointOnLine,direction,inter2)){
        intersections.push_back(inter2[0]);
        intersections.push_back(inter2[1]);
    }
    else{
        return false;
    }
    //Identifico i 4 punti trovati sulla retta
    double a = ((intersections[0]-pointOnLine).dot(direction))/(direction.norm()*direction.norm());
    double b = ((intersections[1]-pointOnLine).dot(direction))/(direction.norm()*direction.norm());
    //a e b sono delle posizioni relative a pointOnLine, ma non so come sono ordinati perché questo dipende
    //dal senso di lettura dei lati dei poligoni e dalla posizione relativa alla retta
    if(abs(a-b)<tol1){ //i due punti di intersezione sono sovrapposti
        return false;
    }
    if(a>b){
        swap(a,b);
        swap(intersections[0],intersections[1]);
    }
    double c = ((intersections[2]-pointOnLine).dot(direction))/(direction.norm()*direction.norm());
    double d = ((intersections[3]-pointOnLine).dot(direction))/(direction.norm()*direction.norm());
    if(abs(d-c)<tol1){
        return false;
    }

    if(c>d){
        swap(c,d);
        swap(intersections[2],intersections[3]);
    }
    if(b-tol1<=c || d-tol1<=a)// i segmenti [a,b] e [c,d] non si intersecano
    {
        return false;
    }
    //c'è un'intersezione "propria", le tracce sono non-passanti per entrambe le fratture
    if(tol1<c-a && tol1<b-c && tol1<d-b)
    {
        E1 = intersections[2];
        E2 = intersections[1];
        tips1 = true;
        tips2 = true;
        return true;
    }
    else if(tol1<a-c && tol1<d-a && tol1<b-d)
    {
        E1 = intersections[0];
        E2 = intersections[3];
        tips1 = true;
        tips2 = true;
        return true;
    }
    //un segmento è contenuto nell'altro, possibilmente con uno o entrambi gli estremi che si sovrappongono
    else if(a-tol1<=c && d<=b+tol1)
    {
        E1 = intersections[2];
        E2 = intersections[3];
        tips2 = false;
        tips1 = true;
        if(abs(a-c)<tol1 && abs(d-b)<tol1)
        {
            tips1 = false;
        }
        return true;
    }
    else if(c-tol1<=a && b<=d+tol1)
    {
        E1 = intersections[0];
        E2 = intersections[1];
        tips1 = false;
        tips2 = true;

        if(abs(a-c)<tol1 && abs(d-b)<tol1)
        {
            tips2 = false;
        }
        return true;
    }
    return false;
}

//***************************************************************************************************************
// Bounding box:
bool near1(const MatrixXd& fracture1, const MatrixXd& fracture2){
    bool check = true;

    Vector3d max1 = fracture1.rowwise().maxCoeff();
    Vector3d min1 = fracture1.rowwise().minCoeff();
    Vector3d max2 = fracture2.rowwise().maxCoeff();
    Vector3d min2 = fracture2.rowwise().minCoeff();

    if (min1(0)> max2(0)+tol1 || min1(1) > max2(1)+tol1 || min1(2) > max2(2)+tol1)
    {
        check = false;
    }
    if (min2(0) > max1(0)+tol1 || min2(1) > max1(1)+tol1 || min2(2) > max1(2)+tol1)
    {
        check = false;
    }

    return check;
}

//***************************************************************************************************************
// Intersezione sfere circoscrtite
bool near2(const MatrixXd& fracture1, const MatrixXd& fracture2){
    bool check = true;
    VectorXd bar_fracture1 = fracture1.colwise().mean();
    VectorXd bar_fracture2 = fracture2.colwise().mean();

    double max_1 = 0.0;
    for (unsigned int j = 0; j < fracture1.cols(); j++){
        double length1 = abs(bar_fracture1.norm() - fracture1.col(j).norm());
        if (length1 > max_1+tol1){
            max_1 = length1;
        }
    }

    double max_2 = 0.0;
    for (unsigned int j = 0; j < fracture2.cols(); j++){
        double length2 = abs(bar_fracture2.norm() - fracture2.col(j).norm());
        if (length2 > max_2+tol1){
            max_2 = length2;
        }
    }

    double control_length = max_1 + max_2;
    double bar_distance = abs((bar_fracture1 - bar_fracture2).norm());

    if (bar_distance + tol1 < control_length){
        cout << "Intersection is possible" << endl;
    }
    else{
        check = false;
        //cout << "Intersection is impossible" << endl;
    }
    return check;
}



//PARTE 2
//***************************************************************************************************************
bool IsInside(const MatrixXd& frattura, const Vector3d& puntoRetta,const Vector3d& direzione,vector<Vector3d>& intersezioni, vector<unsigned int>& lato){
    unsigned int counter = 1;
    unsigned int meet1 = 0;
    unsigned int meet2 = 0;

    for(unsigned int i= 0; i<frattura.cols();i++){
        Vector3d puntoipoly = frattura.col(i);
        Vector3d latoipoly =frattura.col((i+1) % frattura.cols())-frattura.col(i);

        if(!(latoipoly.cross(direzione)).isZero(tol2)){
            double t = ((puntoRetta.cross(direzione)).dot((latoipoly.cross(direzione)))-(puntoipoly.cross(direzione)).dot((latoipoly.cross(direzione))))/(latoipoly.cross(direzione).norm()*latoipoly.cross(direzione).norm());

            if ( (t > (-tol2)) && (t < (1.0+tol2)) ){
                intersezioni.push_back(puntoipoly + t * latoipoly);
                lato.push_back(counter);

                double k = (puntoipoly + t*latoipoly -puntoRetta).dot(direzione)/(direzione.norm()*direzione.norm());
                if(k > tol2){
                    meet1++;
                    if(k<1.0+tol2){
                        meet2++;
                    }
                }
            }
        }
        else{
            if((puntoipoly-puntoRetta).cross(direzione).isZero(tol1)){
                //ATTENZIONE CONTROLLO COMPRESO
                cout << "fakes" << endl;
                lato.push_back(counter);
            }
        }
        counter++; //passo al lato successivo
    }
    if(meet1 == 1 || meet2 == 1){
        return true;
    }
    return false;
}

void split(const MatrixXd& polygon, const vector<unsigned int>& all, const vector<Trace>& tracce, unsigned int counter, list<MatrixXd>& cutted){
    //inizializzazione variabili
    unsigned int current;
    Vector3d point1 = {};
    Vector3d point2 = {};
    Vector3d direzione = {};
    bool full = false;

    vector<Vector3d> intersezioni = {};
    vector<unsigned int> lati;
    lati.reserve(2);

    if(counter < all.size()){
        current = all[counter]; // ID traccia tagliante
        point1 = tracce[current].vertices.col(0); // estremo 1
        point2 = tracce[current].vertices.col(1); // estremo 2
        direzione = point2-point1; // direzione dell traccia

        bool near = near1(tracce[current].vertices,polygon);
        if(near){
            // full = true se la frattura contiene un estremo della traccia
            full = IsInside(polygon,point1,direzione,intersezioni,lati);
            if( lati[0] == lati[1]){
                full = false;
            }
        }

        counter++;
        if( !full){
            split(polygon,all,tracce,counter,cutted);
        }
        else{
            if(intersezioni.size() != 2){
                cout << "error: " << intersezioni.size() << endl;
            }
            else{
                unsigned int dim = polygon.cols();
                MatrixXd temp_sx(3,dim+2);
                MatrixXd temp_dx(3,dim+2);

                bool change = false;
                unsigned int counter1 = 0;
                unsigned int counter2 = 0;

                for(unsigned int i = 0; i< dim ; i++ ){
                    if(!change){
                        temp_sx.col(counter1) = polygon.col(i);
                        counter1 ++;
                    }
                    if(change){
                        temp_dx.col(counter2) = polygon.col(i);
                        counter2 ++;
                    }
                    if( (lati[0]-1) == i){
                        if(intersezioni[0]!=polygon.col(i)){
                            temp_sx.col(counter1) = intersezioni[0];
                            counter1++;
                        }
                        change = true;
                    }
                    if((lati[1]-1) == i){
                        if(intersezioni[1]!=polygon.col(i)){
                            temp_sx.col(counter1) = intersezioni[1];
                            counter1 ++;

                            temp_dx.col(counter2) = intersezioni[1];
                            counter2 ++;
                            temp_dx.col(counter2) = intersezioni[0];
                            counter2 ++;
                        }
                        change = false;
                    }
                }

                // estraggo sotto matrice riempita
                MatrixXd sx = temp_sx.leftCols(counter1);
                MatrixXd dx = temp_dx.leftCols(counter2);

                // chiamate ricorsive
                split(sx, all, tracce, counter, cutted);
                split(dx, all, tracce, counter, cutted);
            }
        }
    }
    else{
        cutted.push_back(polygon);
    }
}

list<MatrixXd> cutting(vector<Fracture>& fratture, vector<Trace>& tracce){ // const
    list<MatrixXd> cutted = {};  // lista con le fratture tagliate (ogni matrice ha i vertici di un taglio)

    // Ciclo su tutte le fratture
    for(unsigned int i = 0; i<fratture.size(); i++){
        MatrixXd polygon = fratture[i].vertices; // vertici frattura

        // Unifico vettore passanti e non passanti
        vector<unsigned int> passing = fratture[i].passing;
        vector<unsigned int> not_passing = fratture[i].not_passing;
        vector<unsigned int> all;
        all.reserve(passing.size() + not_passing.size());
        all.insert(all.end(), passing.begin(), passing.end());
        all.insert(all.end(), not_passing.begin(), not_passing.end());

        // Se la frattura ha delle tracce faccio i tagli
        if(!all.empty()){
            unsigned int counter = 0;
            split(polygon, all , tracce, counter, cutted); // funzione ricorsiva taglio finchè posso
        }
        else{
            cutted.push_back(polygon);
        }
    }

    // // CHECK:
    // for(auto const& i :cutted)
    //     cout << i << endl << endl;

    cout << "tracce: " << tracce.size() << endl;
    cout << "cutted: "<<cutted.size() << endl;

    return cutted;
}

void extractinfo(const list<MatrixXd>& cutted, Mesh& mesh){
    unsigned int dim = cutted.size();

    //inizializzo strutture dati necessarie
    unordered_map<Vector3d, unsigned int, Vector3dHash> map_point;
    vector<unsigned int> C0Id;
    C0Id.reserve(dim*4);
    vector<Vector3d> C0;
    C0.reserve(dim*4);

    unordered_set<pair<unsigned int, unsigned int>,pair_hash> set_edge;
    vector<unsigned int> C1Id;
    C1Id.reserve(dim*4);
    vector<Vector2i> C1V;
    C1Id.reserve(dim*4);

    vector<unsigned int> C2Id;
    C2Id.reserve(dim);
    vector<vector<unsigned int>> C2V;
    C2V.reserve(dim);
    vector<vector<unsigned int>> C2E;
    C2E.reserve(dim);

    unsigned int Idp = 0;
    unsigned int position = 0;

    unsigned int prev = 0;
    unsigned int Ide = 0;
    unsigned int postion_edge = 0;
    unsigned int Ida = 0;

    for (auto const& polygon : cutted) {
        C2Id.push_back(Ida);
        Ida++;
        vector<unsigned int> vertices;
        vector<unsigned int> edges;
        for(unsigned int i = 0; i<polygon.cols()+1; i++){
            Vector3d point = polygon.col(i % polygon.cols());

            int rounding = -static_cast<int>(log10(tol1)); //forza conversione in intero
            double power = pow(10, rounding);
            for (unsigned int i = 0; i < point.size(); ++i) {
                point(i) = round(point(i)*power)/power;
            }

            // Celle0D
            auto result1 = map_point.insert({point,Idp});
            if(result1.second){
                C0Id.push_back(Idp);
                C0.push_back(point);
                position = Idp;
                Idp++;
            }
            else{
                position = result1.first->second;
            }

            //Celle1D
            if(i>=1){
                pair<unsigned int, unsigned int> edge;
                if(prev <= position){
                    edge.first =  prev;
                    edge.second = position;
                }
                else{
                    edge.first = position;
                    edge.second =  prev;
                }

                auto result2 = set_edge.insert(edge);
                if(result2.second){
                    C1Id.push_back(Ide);
                    C1V.push_back({prev,position});
                    postion_edge = Ide;
                    Ide++;
                }
                else{
                    postion_edge =  result2.first->second;
                }

                edges.push_back(postion_edge);
            }

            prev = position;
            if(i<polygon.cols()){
                vertices.push_back(position);
            }
        }
        C2V.push_back(vertices);
        C2E.push_back(edges);
    }

    //CHECk:
    // cout << "punti" << endl;
    // for(unsigned int i = 0; i<C0Id.size(); i++){
    //     cout << C0Id[i] << ": " << C0[i].transpose() << endl;
    // }
    // cout << endl;

    // cout << "lati" << endl;
    // for(unsigned int i = 0; i<C1Id.size(); i++){
    //     cout << C1Id[i] << ": " << C1V[i].transpose() << " (  " << C0[C1V[i][0]].transpose() << ", " << C0[C1V[i][1]].transpose() << " )" << endl;
    // }

    // cout << endl;
    // cout << "aree" << endl;
    // for(unsigned int i = 0; i<C2Id.size(); i++){
    //     cout << C2Id[i] << ": ";
    //     for(unsigned int j = 0; j<C2V[i].size(); j++){
    //         cout << C2V[i][j] << ", ";
    //     }
    //     cout << "    ";
    //     for(unsigned int k = 0; k<C2E[i].size(); k++){
    //         cout << C2E[i][k] << ", ";
    //     }
    //     cout << endl;
    // }

    mesh.NumberCell0D = Idp;
    mesh.Cell0DId = C0Id;
    mesh.Cell0DCoordinates = C0;

    mesh.NumberCell1D = Ide;
    mesh.Cell1DId = C1Id;
    mesh.Cell1DVertices = C1V;

    mesh.NumberCell2D = Ida;
    mesh.Cell2DId = C2Id;
    mesh.Cell2DVertices = C2V;
    mesh.Cell2DEdges = C2E;
}
}




