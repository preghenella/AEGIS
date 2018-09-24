#ifndef ALIGENEPEMV1_H
#define ALIGENEPEMV1_H
/* Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 * Copyright(c) 1997, 1998, 2002, Adrian Alscher and Kai Hencken          *
 * Copyright(c) 2002 Kai Hencken, Yuri Kharlov, Serguei Sadovsky          *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Event generator of single e+e- pair production in ultraperipheral PbPb collisions
// at 5.5 TeV/nucleon.
// Author: Yuri.Kharlov@cern.ch
// 9 October 2002
//
// Revised on September 2018 for ALICEo2: Roberto Preghenella (preghenella@bo.infn.it)

#include "TEpEmGen.h"
#include "TRandom.h"

//-------------------------------------------------------------
class TGenEpEmv1 : public TEpEmGen {
  
 public:
  TGenEpEmv1();
  virtual ~TGenEpEmv1();

  virtual void GenerateEvent();
  virtual void Init();
  void SetDebug(Int_t debug) {fDebug=debug;}
  void SetYRange(Double_t min, Double_t max) {fYMin = min; fYMax = max;};
  void SetPtRange(Double_t min, Double_t max) {fPtMin = min; fPtMax = max;};
  void SetPhiRange(Double_t min, Double_t max) {fPhiMin = min; fPhiMax = max;};
  void SetOrigin(Float_t ox, Float_t oy, Float_t oz) {fOrigin[0]=ox; fOrigin[1]=oy; fOrigin[2]=oz;};
  void SetSigma(Float_t sx, Float_t sy, Float_t sz) {fOsigma[0]=sx; fOsigma[1]=sy; fOsigma[2]=sz;};
  void SetTimeOrigin(Float_t timeorig) {fTimeOrigin = timeorig;};
  
 protected:
  TGenEpEmv1(const TGenEpEmv1 & gen);
  TGenEpEmv1 & operator=(const TGenEpEmv1 & gen);
  void GeneratePair(Double_t vx, Double_t vy, Double_t vz, Double_t vt);
  void Rndm(Float_t *array, Int_t n) {gRandom->RndmArray(n, array);};
  void Rndm(Double_t *array, Int_t n) {gRandom->RndmArray(n, array);};

  Float_t    fMass;    // electron mass
  Int_t      fDebug;   // debug level
  Int_t      fEvent;   // internal event number

  Double_t fYMin, fYMax;
  Double_t fPhiMin, fPhiMax;
  Double_t fPtMin, fPtMax;
  Double_t fTimeOrigin, fOrigin[3], fOsigma[3];

  ClassDef(TGenEpEmv1,1) // Generator of single e+e- pair production in PbPb ultra-peripheral collisions
};
#endif
