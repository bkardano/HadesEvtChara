# Event-by-Event Event Plane Resolution Correction

A lightweight utility for applying per-event resolution corrections to flow observables in HADES analyses.  
Instead of dividing by a single averaged **Rn**, each event is corrected with an individual **Rn** value looked up from a 2D matrix indexed by `(N_TOFRPC, Z_FW)`.

Reference: https://inspirehep.net/literature/1207642

---

## Files

| File | Purpose |
|---|---|
| `HFlowEPres.h` | Class that loads the correction matrix and serves Rₙ values per event |
| `pseudo-code-e-by-e-EP-corr.C` | Reference macro showing how to integrate the correction into an analysis loop |

---

## How It Works

The event plane resolution **Rₙ** depends on two event-level quantities:

- **`N_TOFRPC`** — hit multiplicity in  TOF-RPC
- **`Z_FW`** — total charge of spectator fragments in the Forward Wall

A pre-computed ROOT file contains one `TH2` per harmonic `n`, where each bin holds the resolution factor **Rn** for that `(N_TOFRPC, Z_FW)` combination.

---

## Step-by-Step Usage

### 1. Include the header

```cpp
#include "HFlowEPres.h"
```

### 2. Point to the correction file

Set the path to the `.root` file that contains the 2D resolution histograms.  
The file must exist for the beamtime you are analysing, for example:

for Apr12 Au+Au 1.23 AGeV

```cpp
TString EPresCorrFile = "EPHistQA_apr12_gen8_pass37_EPResolution_TOFRPC_FWCharge_ParameterFileCal.root";
```
or for Mar19 Ag+Ag 1.58 AGeV 

```cpp
EPresCorrFile = "EPHistQA_mar19_ag158ag_3200A_gen6_pass11_EPResolution_TOFRPC_FWCharge_ParameterFileCal.root";
```

or for Mar19 Ag+Ag 1.23 AGeV 

```cpp
EPresCorrFile = "EPHistQA_mar19_ag123ag_2500A_gen5_pass8_EPResolution_TOFRPC_FWCharge_ParameterFileCal.root";
```
      

### 3. Initialise the correction object

```cpp
Bool_t useEPcorr = kTRUE;
Int_t  nHarmonics = 6;

HFlowEPres* EPresCorr = nullptr;
vector<Float_t> EPres(nHarmonics, 1.0);   // default: no correction

if (useEPcorr && EPresCorrFile != "") {
    EPresCorr = new HFlowEPres("HFlowEPres", "HFlowEPres");
    TFile* fileEPresCorr = TFile::Open(EPresCorrFile.Data(), "READ");

    if (fileEPresCorr && fileEPresCorr->IsOpen()) {
        EPresCorr->LoadHistograms(fileEPresCorr);
    } else {
        cerr << "ERROR: Could not open EP resolution file." << endl;
        return 1;
    }
}
```

> `LoadHistograms()` reads all 8 harmonics from the file and stores them internally.  
> It prints a warning and returns `kFALSE` if any histogram is missing.

### 4. In the event loop — retrieve event-level quantities

```cpp
// Event plane angles
Float_t ep  = evtChara->getEventPlane(eEPcorr);
Float_t epA = evtChara->getEventPlane(eEPcorr, 1);   // sub-event A
Float_t epB = evtChara->getEventPlane(eEPcorr, 2);   // sub-event B
if (ep == -1 || epA == -1 || epB == -1) continue;    // skip if EP not reconstructed

// Quantities that index the correction matrix
Int_t TOFRPC      = evtChara->getCentralityEstimator(HParticleEvtChara::kTOFRPC);
Int_t FWSumChargeZ = evtChara->getCentralityEstimator(HParticleEvtChara::kFWSumChargeZ);
```

### 5. Look up Rₙ for this event

```cpp
// Reset to 1 (= no correction) for each event
fill(EPres.begin(), EPres.end(), 1.0);

if (EPresCorr) {
    for (Int_t iHarm = 0; iHarm < (Int_t)EPres.size(); iHarm++) {
        Double_t Rn = EPresCorr->GetEPresCorr(iHarm + 1, TOFRPC, FWSumChargeZ);
        if (Rn > 0.0 && Rn < 0.99) EPres[iHarm] = Rn;
    }
}
```

> `GetEPresCorr(harmonic, N_TOFRPC, Z_FW)` — harmonic index starts at **1**.  
> Returns `-1` if the histogram is missing or the index is out of range.  
> Values outside `(0, 0.99)` are ignored and the default of `1` is kept.

### 6. In the track loop — apply the correction

```cpp
for (Int_t iHarm = 0; iHarm < nHarmonics; iHarm++) {
    Float_t deltaPhi = phi - ep;   // azimuthal angle relative to event plane

    Float_t Vn = TMath::Cos((iHarm + 1.) * deltaPhi);
    Float_t Sn = TMath::Sin((iHarm + 1.) * deltaPhi);

    // Divide by the event-specific resolution factor
    if (EPres[iHarm] > 0.0 && EPres[iHarm] < 1.0) {
        Vn /= EPres[iHarm];
        Sn /= EPres[iHarm];
    }

    // Fill flow profiles
    fProfile_PtYCent_FlowCos[iHarm]->Fill(Y, pt, cent, Vn, weight);
    fProfile_PtYCent_FlowSin[iHarm]->Fill(Y, pt, cent, Sn, weight);
}
```

---

## Validation

Both the event-by-event method and the conventional approach (dividing by a single mean **Rₙ**) yield consistent results for `v₁(pT)` — shapes and absolute values agree within statistical uncertainties.

---


