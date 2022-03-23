// SCDWrapper -- SCDWrapper that manipulate the whole bulk
#pragma once
#include"Damage.h"
#include"Bundle.h"
#include"cpdf.h"
#include"rvgs.h"
#include"CascadeDamage.h"
//#include"constant.h"
//#include"gnuplot_i.h"
#include<string>

class SCDWrapper {
private:
    /* private data member */
    unordered_map<int64, Object*> allObjects;  // map that store all object
    unordered_map<int64, Object*> mobileObjects;  // map that store mobile object
    unordered_map<Object*, Bundle*> linePool;
    unordered_map<int64, int> surface;
    unordered_map<int64, int> bottom;
    Damage damage;
    Cpdf cpdf;
    int sinks[LEVELS+1][POINTS]; /* sinks[0][n] vacancy, sinks[1][n] interstitial*/
    long double sinkDissRate[2][POINTS];
    // dissociation rate of V/H from dislocations
    int reactions[8][POINTS];
    double dislocationD; // [cm^-2] dislocation density
    int round;
    /* hold reactions, 1st dimension is reaction type, second dimension is element, value is total number of this reaction */
    long double matrixRate[POINTS];    // total rate in every element(point)
    long double bulkRate;  // total rate in the whole bulk;
    enum InsertStyle {INTERSTITIAL, SUBSTITUTIONAL};
    ofstream fs1, fs2, fs3, fs4, fs5, fs6;
    fstream fs;
    fstream fv;
    //GnuplotS gs, gr; /* plot species.out and reaction */
    //GnuplotS gd1, gd2; /* damage graph 1 and damage graph 2*/
    //GnuplotS gh1, gh2, gh3; /* H deposition graph 1,2,3 */
    //GnuplotS gv;
    
    /* SCD preparations */
    int iStep;
    double advTime;
    double accTime;
    double addTime;
    double totalDPA;
    double dpa;
    fstream st;
    double one_Over_L;
    double averageR;
    double delta_t;
    double absM;
    double absRho;
    
    /* private functions */
    /* set sinks function */
    void setSinks();
    /* function to compute dissociation reaction rate */
    void computeSinkDissRate(const int, const int);
    /* a way to get key from attribute */
    int64 attrToKey(const int* const); /* same with Object::setKey() */
    /* decide insertion mode function */
    int64 atomProperty(InsertStyle, const int&);
    /* moniter map functions */
    void addNewObjectToMap(const int64&, const int*);
    /* add to map by key(increase m to nth element)---> this is created for restart*/
    void addNewObjectToMap(const int64&, const int&);
    /* add new object to map by key (only increase by 1 to #kth element)*/
    void addNewObjectToMap(Object*);
    /* add new object to map by object pointer(only increase 1 to #nth element) */
    void updateObjectInMap(Object*, const int&);
    void removeObjectFromMap(const int64&); /* remove one object from map */
    void addReactionToOther(const Object* const);
    /* Impact of one new mobile object to other existing objects */
    void updateRateToOther(const Object* const, const int&);
    /* when number of this object changes, rates related to this object change also */
    void removeRateToOther(const int64&);
    /* when one mobile object has been removed,
     ** rates of this object with other extisting objects should also be removed
     */
    void updateSinks(const int, const int*); /* only for restart use */
    /* process event functions */
    void processDiffEvent(Object*, const int&, const char&);     /* process diffusion reactionObject */
    void processSinkEvent(Object*, const int&);     /* process absorption reaction */
    void processDissoEvent(Object*, const int&, const int64&, fstream& ); /* process dissociation event */
    void processCombEvent(Object*, const int&, const int64&, fstream& );  /* process combination reaction */
    void processSinkDissEvent(const int, const int&); /* process dissociation from sink event */
    /* get insertion functions */
    void getElectronInsertion(const int&);
    void getNeutronInsertion(const int&);
    void getIonInsertion(const int&, const double&, fstream&);
    void getParticleInsertion(const int&, const double&, fstream&);    // deal with damage[0]
    void getHeInsertion(const int&);  // deal with damage[1]
    void getHInsertion(const int&, const double&, fstream&);   // deal with damage[2]
    /* write file funcitons*/
    void writeSinkFile(const Object* const, const int& n, const double&);
    /* sinks.out only writes when sink reaction happens, now this function is not "writing things" but only updating sinks[][] */
    void writeSpeciesFile(const double&, const int&);
    void writeClusterFile(const double&, const int&);
    /* distinguish reaction type function */
    bool recognizeSAV(const Object* const, const Object* const);
    /* if SAV, return true, if not, return false*/
    int countDefectNumber(const int &, char*);
    /* this function counts the number of defects(V,SIA,H...) in every element, returns total number of this defect in the bulk */
    void countRatioDistribution(double&);
    /* count ratio of H to V in the bulk and at the same time decide whether to supplement H to system */
    void sizeDistribution(double&); /* get size distribution */
    void writeReaction(); /* take down reactions for drawing */
    
public:
    SCDWrapper();  // constructor: for start ;
    //SCDWrapper();  // constructor: for restart;
    void computeMatrixRate(const int& n);  // computes total rate in element n
    void computeBulkRate();
    Object* selectReaction(int64&, Reaction&, int&);  // select reaction
    void processEvent(const Reaction&, Object*, const int&, const int64&, const double&, const double&);    // deal with reactions
    ~SCDWrapper();          /* destructor to delete everything newed */
    // get series functions that allow direct manipulation on private data member
    unordered_map<int64, Object*>* getAllObjects();
    unordered_map<int64, Object*>* getMobileObjects();
    unordered_map<Object*, Bundle*>* getLinePool();
    const double getAndExamineRate(); /* computes bulk rate and matrix rate and returns bulk rate*/
    /* output file functions */
    void displayDamage();
    void displayAllObject();
    void writeFile(const double&, const int&);
    /* testing functions which includes display/output file functions that might be used when debugging */
    void writeSinkFile(const double&, const int&);
    void display() const;   // display rateMatrix within in one element
    void restart(int&, double&);
    /* test functions, don't need in the future */
    void test(const int&); /* count how many vacancies now in the bulk */
    
    /* draw picture functions */
    void drawSpeciesAndReactions(double&t);
    void drawDamage(double&); /* this function draw all diagrams in damage process */
    void drawHD(double&); /* this function draw all diagrams in H deposition process */
    void writeVacancy();
    double getTotalDpa(double&);
    void doClusterDynamics(double&, double&, double&, double&, double&, double&); // this function does one round of SCD and transfers values of oneOverL and dt to CP
    void updateRates(double&); // CP send in renewed dislocation density, and the SCD should update all rates according to the renewed information.
    void clearSRSCD();
};

