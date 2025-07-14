# FinFlutterHandbook
A guide for predicting the critical flutter speed of a rocket fin.

Aeroelastic theory taken from [John K. Bennett's article in Peak of Flight issue 615](https://www.apogeerockets.com/Peak-of-Flight/Newsletter615), which corrected errors present in Peak of Flight issues 291 and 411. More details on his work can be found in his repository [here](https://github.com/jkb-git/Fin-Flutter-Velocity-Calculator/tree/main). This work puts his fin flutter calculations into Python functions and combines them with an atmospheric model that better represents the conditions of specific launch sites.

This repository contains:
- A Python function for determining the critical flutter speed of a trapezoidal fin according to Bennett's corrected theory: `flutter_velocity_trapezoidal_fin`
- Code for determining the critical flutter speed of an arbitrary polygonal fin according to Bennett's corrected theory (in progress)
- An explanation of the theory behind fin flutter and other methods for determining the critical flutter speed  (in progress)

To be added:
- Excel sheet with the calculations and ability to define location-specific atmospheric conditions
- Code for determining the critical flutter speed of an elliptical fin


links to open next time (not sure all useful):
https://sdasoftware.com/software/nastran/features/aeroelasticity/
https://www.nakka-rocketry.net/RD_fin.htm
https://rocketry.gitbook.io/public/tutorials/airframe/sizing-fins
https://www.mas.bg.ac.rs/_media/istrazivanje/fme/vol50/3/17_m._dinulovic_et_al.pdf
https://www.ata-e.com/services/analysis/aeroelasticity-flutter/
https://software.nasa.gov/software/DRC-014-036
https://www.apogeerockets.com/Peak-of-Flight/Newsletter615
https://www.reddit.com/r/AerospaceEngineering/comments/q9ml9f/aeroelastic_flutter_analysis/


check that these are using outdated models:
https://github.com/wisconsinspaceprogram/aero-analysis/tree/main
https://github.com/steelmetro2055/Fin-Flutter/tree/main
https://github.com/wisconsinspaceprogram/aero-analysis
    nice for ork integration? Do something similar with this repo?
https://github.com/UTATRocketry/Flutter-Analysis