#include "cpdf.h"

Cpdf::Cpdf()
{
    
    char cpdf[20];
    ifstream fc;
    /*initialize max Possibility*/
    for(int i=0; i<POINTS; i++){
        maxPossibility[i] = 0;
    }
    for(int fileNumber = 0; fileNumber < POINTS; fileNumber++){
    	/* first point is surface */
        vector<double> tempEnergy, tempCumul;
        double energy= 0.0 , cumul = 0.0;
        sprintf(cpdf,"cpdf%d.txt",fileNumber);
        fc.open(cpdf);
        string oneLine;
        stringstream lineHold;
        while (getline(fc, oneLine)) {
            lineHold.str(oneLine);
            lineHold >> energy >> cumul;
            tempEnergy.push_back(energy);
            tempCumul.push_back(cumul);
            lineHold.clear();
        } //read a file and store/create corresponding vectors
        pair<int, vector<double> > newNodeOne(fileNumber, tempEnergy);
        pair<int, vector<double> > newNodeTwo(fileNumber, tempCumul);
        recoilEnergy.insert(newNodeOne);
        cumulPossibility.insert(newNodeTwo);
        /* add new nodes to members */
        size[fileNumber] = (int)tempEnergy.size();
        maxPossibility[fileNumber] = tempCumul[(int)tempEnergy.size()-1];
        /* store size and possibility for SCD usage in case */
        tempEnergy.clear();
        tempCumul.clear();
        /* clean two vectors and start next loop */
        fc.close();
    }
}

double Cpdf::samplePkaEnergy(const double & xi, const int & n) const
{
    int  i, index = 0;
    double slope;
    unordered_map<int, vector<double> >::const_iterator gotEnergy = recoilEnergy.find(n);
    unordered_map<int, vector<double> >::const_iterator gotCumul = cumulPossibility.find(n);
    vector<double> energy = gotEnergy->second;
    vector<double> cumul = gotCumul->second;
    for (i = 1; i < size[n] ; ++i) {
        if ((xi > cumul[i - 1]) && (xi < cumul[i])) {
            index = i;
            break;
        }
    }
    /* Interpolation: */
    slope = (energy[index] - energy[index - 1]) / (cumul[index] - cumul[index - 1]);
    return slope*(xi - cumul[index - 1]) + energy[index - 1];
}

double Cpdf::getMaxPossibility(const int& n){
    return maxPossibility[n];
}
