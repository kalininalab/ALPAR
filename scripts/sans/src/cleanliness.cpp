#include "cleanliness.h"


void cleanliness::init() {
}

void cleanliness::addWeightStateBefore(double weight, color_t color) {
    struct filter_weight_stats st = {};
    st.weight = weight;
    st.color = color;
    this->weightsBefore.push_back(st);
}

void cleanliness::setSmallestWeight(double weight, color_t color) {
    this->smallestWeight = {};
    this->smallestWeight.weight = weight;
    this->smallestWeight.color = color;
    this->incrementWeightAfterCounter(weight);

}

void cleanliness::calculateWeightBeforeCounter() {
    bool smallestReached = false;

    for (int i = 0; i< weightsBefore.size() && !smallestReached; i++) {
        this->splitCountBefore++;
        this->splitWeightCountBefore+= this->weightsBefore[i].weight;

        if (this->weightsBefore[i].is(this->smallestWeight)) {
            smallestReached = true;
        }
    }
}

void cleanliness::incrementWeightAfterCounter(double weight) {
    this->splitWeightCountAfter+=weight;
}

void cleanliness::setFilteredCount(uint64_t splitCountAfter) {
    this->splitCountAfter = splitCountAfter;
}

void cleanliness::reportCleanliness() {
    cout << "cleanliness: " <<  this->splitCountAfter << "/" <<  this->splitCountBefore;
    if (this->splitCountBefore > 0) {
        cout  << " -> " << (this->splitCountAfter /  (double)  this->splitCountBefore) << endl;
    } else {
        cout << endl;
    }

    cout << "weighted cleanliness: " <<  this->splitWeightCountAfter << "/" <<  this->splitWeightCountBefore;
    if (this->splitWeightCountBefore > 0) {
        cout << " -> " << (this->splitWeightCountAfter / this->splitWeightCountBefore) << endl;
    } else {
        cout  << endl;
    }
}