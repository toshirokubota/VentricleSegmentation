#ifndef ___SIZE_FEATURE_DEFAULT_PARAM_H___
#define ___SIZE_FEATURE_DEFAULT_PARAM_H___

const int SZ_LABEL_NO_NODULE = 0;
const int SZ_LABEL_SOLID_NODULE = 1;
const int SZ_LABEL_NONSOLID_NODULE = 2;
const int SZ_LABEL_PARTSOLID_NODULE = 3;
const int SZ_LABEL_NON_BACKGROUND = 4;
const int SZ_LABEL_SOLID_NODULE_LOW = 10001; //fuzzy solid due to canonical sub-sampling
const int SZ_LABEL_NONSOLID_NODULE_LOW = 10002; //fuzzy non-solid due to canonical sub-sampling
const int SZ_LABEL_SOLITARY_NODULE = 10;
const int SZ_LABEL_JUXTAPOSITION_NODULE = 20;

const int DefaultMaxSolitarySize = 500;
const int DefaultRDNumIteration = 4;
const float DefaultRDLambda = (float)0.5;
//const float DefaultFindSeedThres = (float)10.0;
const int DefaultFindSeedWidth = 3;
const int DefaultFindSeedCoreWidth = 3;
const float DefaultFindSeedPerc = (float)0.95;
const int DefaultMaxNumSeeds = 1;
const float DefaultFindSeedThreshold = 0.5; //allows two neighbors to be the same value
//const float DefaultFindSeedThresholdOld = 25.0; //allows two neighbors to be the same value

const float DefaultMinimumValidSphericity = 0.1; //anything below is ignored in seed selection.
const float DefaultTraceDownThres = 0;

const int ForegroundColor = 255;

const float DefaultCanonicalVoxelSize[] = {0.85, 0.85, 0.85};

#endif /* ___SIZE_FEATURE_DEFAULT_PARAM_H___ */