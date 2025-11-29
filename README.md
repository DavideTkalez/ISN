# EMG Signal Acquisition and Analysis – Neuromuscular System Engineering (Politecnico di Torino)

This repository contains the material developed for the project of the course
**Neuromuscular System Engineering**, dedicated to EMG signal acquisition and
processing through two experimental protocols: one focused on myoelectric
classification and one on muscle characterization using electrode arrays.

---

## 📁 Repository Contents
- MATLAB scripts for signal preprocessing, feature extraction, classification, and spectral analysis.
- Generated plots (PSD, conduction velocity, fatigue analysis).
- Documentation of the processing workflow.

---

## 🧪 Protocol 1 – Myoelectric Classification (Forearm)

EMG signals were recorded from the forearm using circular electrodes while the subject
performed:
- Hand opening  
- Hand closing

### Goals
- Filter raw EMG and extract envelopes
- Build m × n datasets (samples × channels)
- Assign class labels and split data (70–30)
- Train and evaluate:
  - **Support Vector Machine (SVM)**
  - **Linear Discriminant Analysis (LDA)**
  - **Cosine-similarity classifier** based on ARV prototypes
- Compare classification performance across methods

---

## 🧪 Protocol 2 – Biceps/Triceps Analysis with Electrode Arrays

High-density EMG was collected using two 8-electrode arrays during isometric contractions
with increasing weight loads (e.g., 2–8 kg).

### Goals
- Compute monopolar, single-differential (SD), and double-differential (DD) signals
- Estimate **Power Spectral Density (PSD)** with 250 ms epochs
- Compute **conduction velocity (CV)** using different methods
- Generate fatigue plots (ARV, RMS, MDF, MNF, CV)
- Create linear trend plots across contractions
- Simulate **crosstalk** by summing biceps and triceps recordings
- Recompute PSD and CV to evaluate crosstalk effects
- Apply **source separation** methods to mitigate crosstalk

---

## 🛠 Tools & Technologies
- **MATLAB** for signal processing and classification
- EMG16 and Cometa amplifiers (PoliTO laboratory equipment)

---

## 👤 Authors
Project completed for the course *Neuromuscular System Engineering*  
Politecnico di Torino – Academic Year 2023/2024  
Supervisors: Prof. **Marco Mesin**, Prof. **Matteo Raggi**, and Prof. **Taian Martins Vieira**.

---

## 📄 License
Specify the license here (MIT, GPL, etc.).
