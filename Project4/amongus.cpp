#include <iostream>
#include <string>
#include <sstream>
#include <cmath> 
#include <set>
#include <stack> 
#include <utility>
#include <deque>
#include <queue>
#include <vector>
#include <getopt.h>
#include <fstream>
#include <stdint.h>
#include <cmath>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <limits>
#include <iomanip>
struct Coordinates {

public:

    enum class Location {Lab, Decon, OpenArea};

    int xCoord;
    int yCoord;
    char loc;

    Coordinates(int x_in, int y_in, char loc_in)
        : xCoord(x_in), yCoord(y_in), loc(loc_in) {}

};

struct Best {
public: 
    std::vector<uint32_t> order;
    double value;
};

struct Primfo {
public:
    bool k;
    double d;
    uint32_t p;

    Primfo(bool k_in, double d_in, uint32_t p_in)
        : k(k_in), d(d_in), p(p_in) {}

    Primfo()
        : k(false), d(std::numeric_limits<double>::infinity()), p(0) {}
};



void findDistance(double& curDistance, size_t cord1, size_t cord2, std::vector<Coordinates>& stuff) {
    if ((stuff[cord1].loc == 'l' && stuff[cord2].loc == 'o') || (stuff[cord1].loc == 'o' && stuff[cord2].loc == 'l')) {
        curDistance = std::numeric_limits<double>::infinity();
        return;
    }
    else {
        double yDif = (double)stuff[cord1].yCoord - (double)stuff[cord2].yCoord;
        yDif *= yDif;

        double xDif = (double)stuff[cord1].xCoord - (double)stuff[cord2].xCoord;
        xDif *= xDif;

        curDistance = sqrt(yDif + xDif);
        return;
    }


}


double ghostDistance(Coordinates first, Coordinates second) {
    double xDif = first.xCoord - second.xCoord;
    double yDif = first.yCoord - second.yCoord;

    return sqrt(xDif * xDif + yDif * yDif);
}

std::string mode;

void getMode(int argc, char* argv[]) {
    // bool modeSpecified = false;

    opterr = false;
    int choice;
    int option_index = 0;
    option long_options[] = {
        { "mode", required_argument, 0, 'm'},
        { "help",  no_argument,  nullptr, 'h'},
        { nullptr, 0, nullptr, '\0'}
    };
    //  Fill in the double quotes, to match the mode and help options.
    while ((choice = getopt_long(argc, argv, "m:h",
        long_options, &option_index)) != -1) {
        switch (choice) {
        case 'm':
            mode = optarg;
            break;

        case 'h':
            std::cout << "This program is used to help Marco find the Countess lost in the caslt\n";
            exit(1);
            break;

        default:
            std::cout << "Unkown command line option";
            exit(1);
        }
    }

}

void initVector(std::vector<Coordinates>& points) {
    size_t numPoints;
    std::cin >> numPoints;
    int x, y;
    char loc;

    for (size_t i = 0; i < numPoints; ++i) {
        std::cin >> x >> y;
        if (x < 0 && y < 0)
            loc = 'l';
        else if ((x == 0 && y < 0) || (y == 0 && x < 0) || (y == 0 && x == 0))
            loc = 'd';
        else
            loc = 'o';

        Coordinates curCoords(x, y, loc);
        points.push_back(curCoords);
    }
}

void mst(std::vector<Coordinates> points) {
    using namespace std;
    double totalWeight = 0;

    vector<Primfo> primTable;

    //Initialize vector, could be better way to do this
    for (size_t i = 0; i < points.size(); ++i) {
        Primfo cur;
        primTable.push_back(cur);
    }

    primTable[0].d = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        uint32_t shortestIndex = 0;
        double smallestDistance = std::numeric_limits<double>::infinity();

        //From all FALSE values, pick smallest distance
        for (uint32_t j = 0; j < points.size(); ++j) {
            if (!primTable[j].k && primTable[j].d < smallestDistance) {
                shortestIndex = j;
                smallestDistance = primTable[j].d;
            }

        }
        //Set to true
        primTable[shortestIndex].k = true;
        totalWeight += primTable[shortestIndex].d;

        double distance;
        
        //For each connecting verticy, change distance accordingly
        for (size_t j = 0; j < points.size(); ++j) {
            if (!primTable[j].k) {
                findDistance(distance, shortestIndex, j, points);
                if (distance < primTable[j].d) {
                    primTable[j].d = distance;
                    primTable[j].p = shortestIndex;
                }
                
            }
        }
    }

    if (totalWeight == numeric_limits<double>::infinity()) {
        cerr << "Cannot construct MST";
        exit(1);
    }

    cout << totalWeight << "\n";
    for (uint32_t j = 1; j < primTable.size(); ++j) {
        cout << min(j, primTable[j].p) << " " << max(j, primTable[j].p) << "\n";
    }

}

uint32_t findClosest(Coordinates curPoint, std::vector<Coordinates> list) {
    uint32_t ind =0;

    double closest = std::numeric_limits<double>::infinity();

    for (uint32_t i = 0; i < list.size(); ++i) {
        double distance = ghostDistance(curPoint, list[i]);

        if (distance < closest && distance != 0) {
            closest = distance;
            ind = i;
        }
    }

    return ind;
}

uint32_t findClosest(double &weight, Coordinates curPoint, std::vector<Coordinates> list) {
    uint32_t ind = 0;

    double closest = std::numeric_limits<double>::infinity();

    for (uint32_t i = 0; i < list.size(); ++i) {
        double distance = ghostDistance(curPoint, list[i]);

        if (distance < closest && distance != 0) {
            closest = distance;
            ind = i;
        }
    }

    weight += closest;
    return ind;
}


Best fastTSP(std::vector<Coordinates> points) {
    using namespace std;

    vector<uint32_t> path;
    Best test;

    //path.resize(points.size());
    double weight = 0;
    path.push_back(0);

    uint32_t nearestInd = findClosest(weight, points[0], points);
    path.push_back(nearestInd);

    //Initialization - Find next vertex i for which Cost of i, j is minimum
    for (uint32_t i = 1; i < points.size(); ++i) {
        //3 steps to Nearest Insertion

        
        double smallest = numeric_limits<double>::infinity();

        for (uint32_t j = 0; j < path.size() - 1; ++j) {

            double d13 = ghostDistance(points[path[j]], points[i]);
            double d23 = ghostDistance(points[path[j + 1]], points[i]);
            double d12 = ghostDistance(points[path[j]], points[path[j + 1]]);

            double curDistance = d13 + d23 - d12;
            if (curDistance < smallest) {
                smallest = curDistance;
                nearestInd = j + 1;
            }
        }

        weight += smallest;
        


        path.insert(nearestInd + path.begin(), i);
        //Selection arbitrarily select a city that is not in loop yet

        //Insertion find edge {i, j} in the loop that minimizes Cik + Ckj - Cij
        
        //insert k between i and j
    }

    weight += ghostDistance(points[path[0]], points[path[path.size()-1]]);
    test.value = weight;
    test.order = path;

    if (mode != "OPTTSP") {
        cout << weight << "\n";
        for (uint32_t i = 0; i < path.size() - 1; ++i) {
            cout << path[i] << " ";
        }
        cout << "\n";
    }
    return test;
}


class Optimal {
public:

    bool inFirstNElements(uint32_t el, size_t n) {
        for (uint32_t i = 0; i < n; ++i) {
            if (curOrder[i] == el)
                return true;
        }
        return false;
    }

    
    void genPerms(size_t permLength) {
        if (permLength == curOrder.size()) {
            // Do something with the path
            curWeight += ghostDistance(points[curOrder[0]], points[curOrder[curOrder.size() - 1]]);
            if (curWeight <= bestWeight) {
                bestWeight = curWeight;
                bestOrder = curOrder;
            }
            curWeight -= ghostDistance(points[curOrder[0]], points[curOrder[curOrder.size() - 1]]);
            return;
        }  // if ..complete path

        if (!promising(permLength)) {
            return;
        }  // if ..not promising

        for (size_t i = permLength; i < curOrder.size(); ++i) {
            std::swap(curOrder[permLength], curOrder[i]);
            //curWeight += ghostDistance(points[curOrder[permLength]], points[curOrder[permLength - 1]]);
            curWeight += distanceMatrix[curOrder[permLength]][curOrder[permLength - 1]];

            genPerms(permLength + 1);
            //for (int j = 0; j < curOrder.size(); ++j) {
             //   std::cout << curOrder[j] << " ";
           // }
            //std::cout << std::endl;
            curWeight -= distanceMatrix[curOrder[permLength]][curOrder[permLength - 1]];

            //curWeight -= ghostDistance(points[curOrder[permLength]], points[curOrder[permLength - 1]]);

            std::swap(curOrder[permLength], curOrder[i]);
        }  // for ..unpermuted elements
    }  // genPerms()

    double prim(size_t permLength) {
        using namespace std;
        double weight = 0;
        int indToStart = 0;
        vector<Primfo> primTable;

        //Initialize vector, could be better way to do this
        for (uint32_t i = 0; i < points.size(); ++i) {
            Primfo cur;

            //SHOULD make it so every vector that's already in the path isn't in the MST
            if (inFirstNElements(i, permLength)) {
                cur.k = true;
                cur.p = 99;
            }
            else
                indToStart = i;

            primTable.push_back(cur);
            //indToStart = i;
        }
        
        primTable[indToStart].d = 0;
        //cout << "Perm length is " << permLength << endl;
        for (size_t i = 0; i < points.size(); ++i) {
            uint32_t shortestIndex = indToStart;
            double smallestDistance = std::numeric_limits<double>::infinity();

            //From all FALSE values, pick smallest distance
            for (uint32_t j = 0; j < points.size(); ++j) {
                if (!primTable[j].k && primTable[j].d < smallestDistance) {
                    shortestIndex = j;
                    smallestDistance = primTable[j].d;
                }

            }
            //Set to true
            primTable[shortestIndex].k = true;
            weight += primTable[shortestIndex].d;

            double distance = numeric_limits<double>::infinity();

            //For each connecting verticy, change distance accordingly
            for (size_t j = 0; j < points.size(); ++j) {
                if (!primTable[j].k) {
                    distance = distanceMatrix[shortestIndex][j];
                    //distance = ghostDistance(points[shortestIndex], points[curOrder[permLength+j]]);
                    //cout << "\nDistance between point " << shortestIndex << " and " << j << " was " << distance << endl;
                    //findDistance(distance, shortestIndex, j, points);
                    if (distance < primTable[j].d && distance != 0) {
                        primTable[j].d = distance;
                        primTable[j].p = shortestIndex;
                    }

                }
            }
        }
        //std::cout << "Weight of primming was " << weight;
        return weight;
    }

    Optimal(std::vector<Coordinates>& stuff) {
        points = stuff;
        bestWeight = std::numeric_limits<double>::infinity();
        curWeight = 0;
    }

    double hypotheticalMST(size_t permLength) {
        double total = 0;

        //Find which points to put into potential mst
        total += prim(permLength);
        return total;

    }

    double connectingBlue(size_t permLength) {
        double totalWeight = 0;

        //Find connecting arm from First to unvisited
        double minDistance = std::numeric_limits<double>::infinity();
        size_t mstSize = curOrder.size() - permLength;

        for (uint32_t i = 0; i < mstSize; ++i) {
            double curDistance;
            curDistance = distanceMatrix[curOrder[0]][curOrder[permLength + i]];
            //curDistance = ghostDistance(points[curOrder[0]], points[curOrder[permLength + i]]);

            if (curDistance < minDistance) {
                minDistance = curDistance;
            }
        }
        totalWeight += minDistance;

        //Find weight for Second blue arm, last in the path so far
        minDistance = std::numeric_limits<double>::infinity();

        for (uint32_t i = 0; i < mstSize; ++i) {
            double curDistance;
            curDistance = distanceMatrix[curOrder[permLength - 1]][curOrder[permLength + i]];
            //curDistance = ghostDistance(points[curOrder[permLength - 1]], points[curOrder[permLength + i]]);

            if (curDistance < minDistance) {
                minDistance = curDistance;
            }
        }

        totalWeight += minDistance;

        //std::cout << " Weight of blue lines was " << totalWeight << std::endl;
        return totalWeight;
    }

    bool promising(size_t permLength) {
        using namespace std;
        if (curWeight > bestWeight)
            return false;

        if (curOrder.size() - permLength <= 5) {
            return true;
        }

        double potential = curWeight;
        
        /*
        double armLength, mstCost, totalEst;

        bool promise = true;

        mstCost = prim(permLength);
        armLength = connectingBlue(permLength);
        totalEst = armLength + mstCost + curWeight;
        if (totalEst > bestWeight)
            promise = false;
        
        
        for (size_t i = 0; i < curOrder.size(); ++i)
            cerr << setw(2) << curOrder[i] << ' ';
        cerr << setw(4) << permLength << setw(10) << curWeight;
        cerr << setw(10) << armLength << setw(10) << armLength;
        cerr << setw(10) << mstCost << setw(10) << totalEst << "  " << promise << '\n';
        */
        //potential += hypotheticalMST();
        potential += prim(permLength);
        potential += connectingBlue(permLength);

        /*
        
        for (int p = 0; p < curOrder.size(); ++p) {
            std::cout << curOrder[p] << " ";
        }
        std::cout << std::endl;
        std::cout << "The low estimate for the remaining graph is " << potential - curWeight << " Meaning the total weight would be " << potential << "\n";
        std::cout << "The current best weight is " << bestWeight << std::endl;*/
        //return true;
        //return promise;
        if (potential <= bestWeight) {
            return true;
        }
        else {
          //  std::cout << "NOT PROMISING\n";
            return false;
        }
    }
    
    void makeDistance() {
        distanceMatrix.resize(points.size());

        for (uint32_t i = 0; i < points.size(); ++i) {
            for (uint32_t j = 0; j < points.size(); ++j) {
                distanceMatrix[i].push_back(ghostDistance(points[i], points[j]));
            }
        }
    }

    void findOPT() {
        using namespace std;
        //Recursively go through every permutation to find the best potential route
        //utilizes pruning to cut out routes that do not look promising.
        Best init = fastTSP(points);
        
        bestWeight = init.value;
        bestOrder = init.order;
        for (uint32_t i = 0; i < points.size(); ++i) {
            curOrder.push_back(i);
        }

        makeDistance();

        genPerms(1);
        //std::cout << "HUH\n";
        std::cout << bestWeight << "\n";
        for (uint32_t i = 0; i < bestOrder.size(); ++i) {
            std::cout << bestOrder[i] << " ";
        }
    }

private:
    std::vector<Coordinates> points;
    std::vector<uint32_t> curOrder;
    std::vector<std::vector<double>> distanceMatrix;
    //Should be in tsp already
    //std::vector<uint32_t> path;
    //Needs to be added to tsp, could just keep track of these with curOrder vector and perm length lol
    //std::vector<uint32_t> notUsed;
    std::vector<uint32_t> bestOrder;
    double bestWeight;
    double curWeight;
};

void optTSP(std::vector<Coordinates>& things) {
    Optimal opt(things);
    opt.findOPT();
}

int main(int argc, char* argv[]) {
    using namespace std;
    cerr << fixed << showpoint << setprecision(2);
    cout << setprecision(2); //Always show 2 decimal places
    std::cout << std::fixed; //Disable scientific notation for large numbers
    std::ios_base::sync_with_stdio(false);
    getMode(argc, argv);

    std::vector<Coordinates> points;

    initVector(points);
    if (mode == "MST")
        mst(points);
    else if (mode == "FASTTSP")
        fastTSP(points);
    else if (mode == "OPTTSP")
        optTSP(points);
    else
        std::cerr << "Unrecognized mode option: " << mode << "\n";

}