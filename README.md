# MI_SJTU_EEG
This repository now contains the Python/Matlab development during master degree. The original EEG-EMG connectivity project has been transferred to an organization repository to involve more contribution.

## Aspects
This project aims to study human-machine mutual learning mechanism including EEG processing and feature extraction and protocol design. Current Current research focuses on graph theory analysis of EEG-EMG functional connectivity. 

## Structure
- Model selection and evaluation of **Motor Imagery Recognition** is performed in Master_Motor_imaginary
- Multi-channels EEG-EMG analysis platform based on *mne* where functional connectivity functions and visualization methods are proposed in *EEG-EMG_platform*. Notably, in this folder, *EEG-EMG problem shooting* aims to break down the coherence calculation and interpret its implication. *firstCoRegistration_eeg* define the preprocessing pipeline and provides visualization tools.
- *paradigm_design* provides experiment instruction and GUI.
- In the Matlab folder, the .m file is an example script for Fuzzy model reference learning control.
- A cleaning robot simulation platform has been developped in Pygame in the cleaning_robot_simulation folder.
- The folder NN in control theory contains jupyter notebook of implementation of RBF-NN and RNN for a control system. Please refer to the report for more information.
