
# AVATREE: An open-source computational modelling framework modelling Anatomically Valid Airway TREE


Citation: Nousias S, Zacharaki EI, Moustakas K (2020) AVATREE: An open-source computational modelling framework modelling Anatomically Valid Airway TREE conformations. PLoS ONE 15(4): e0230259. https://doi.org/10.1371/journal.pone.0230259

Editor: Fang-Bao Tian, University of New South Wales, AUSTRALIA

Received: September 4, 2019; Accepted: February 25, 2020; Published: April 3, 2020

Copyright: © 2020 Nousias et al. This is an open access article distributed under the terms of the Creative Commons Attribution License, which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.

Data Availability: The data underlying the results presented in the study are available from https://vessel12.grand-challenge.org/. The repository is now publicly available at https://gitlab.com/LungModelling/avatree. Furthermore, the outcomes of the presented pipeline are available at https://www.kaggle.com/vvrlabeceupatras/pone-avatree-results.

Funding: This work has been co‐financed by the European Regional Development Fund of the European Union and Greek national funds through the Operational Program Competitiveness, Entrepreneurship and Innovation, under the call RESEARCH – CREATE - INNOVATE Take-A-Breath, under grant agreement No. T1EDK-03832). The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.






![Pipeline](https://gitlab.com/LungModelling/avatree/-/raw/master/figures/pipelinev5.png) 









## Introduction

This paper presents a modelling approach that produces Anatomically Valid Airway tree conformations (AVATREE). Such conformations are adapted to personalized geometry and boundary conditions derived from diagnostic imaging and well-established airway extraction methods. Specifically, this study aims to provide an open-source simulation framework to (i) exploit imaging data so as to provide patient-specific representations (ii) perform structural analysis (iii) extend the segmented airway tree to predict the airway branching across the whole lung volume (iv) visualize probabilistic confidence maps of generation data (v) simulate bronchoconstriction to (vi) access patient-specific airway functionality (vii) perform fluid dynamics simulation in patient-specific boundary conditions to access pulmonary function.





## Requirements

*   QT
*   Boost
*   Zlib
*   CGAL
*   VTK
*   flann
*   PCL
*   ITK
*   FAST Framework


## Imaging data analysis

Employ FAST framework to extract input data.

FAST (Framework for Heterogeneous Medical Image Computing and Visualization) is an open-source cross-platform framework with the main goal of making it easier to do processing and visualization of medical images on heterogeneous systems (CPU+GPU).

FAST has been described in the following research articles. If you use this framework for research please cite them:

*[FAST: framework for heterogeneous medical image computing and visualization](http://www.eriksmistad.no/wp-content/uploads/FAST_framework_for_heterogeneous_medical_image_computing_and_visualization.pdf)  
Erik Smistad, Mohammadmehdi Bozorgi, Frank Lindseth  
International Journal of Computer Assisted Radiology and Surgery 2015*

*[High Performance Neural Network Inference, Streaming, and Visualization of Medical Images Using FAST](https://www.eriksmistad.no/wp-content/uploads/High-Performance-Neural-Network-Inference-Streaming-and-Visualization-of-Medical-Images-Using-FAST.pdf)  
Erik Smistad, Andreas Østvik, André Pedersen  
IEEE Access 2019*


![snapshotCT0001](https://gitlab.com/LungModelling/avatree/-/raw/master/figures/snapshotCT0001.png) 
![snapshotCT0002](https://gitlab.com/LungModelling/avatree/-/raw/master/figures/snapshotCT0002.png) 
![snapshotCT0003](https://gitlab.com/LungModelling/avatree/-/raw/master/figures/snapshotCT0003.png) 


```
airwaySegmentation.exe --filename VESSEL12_08.mhd --export_centerline VESSEL12_08_Centerline.vtk --export_segmentation VESSEL12_08_Segmentation.mhd --seed x,y,z
```

![extraction-of-airway-surface-n-centerline-2](https://gitlab.com/LungModelling/avatree/-/raw/master/figures/extraction-of-airway-surface-n-centerline-2.png) 


```
lungSegmentation.exe  --input VESSEL12_08.mhd --output VESSEL12_08_Lungs.mhd
```


## Generation of extended 1-D representation




```
    bool buildVolumes = true;
	bool buildCenterline = true;
	bool build1DModel = true;
	bool surfaceSampling = false;
	bool buildNormals = false;
	bool build3DModel = false;
	bool refinements = false;
	bool segmentation = false;
	int depth = 12;
	int volumeDepth = 16;
	int density = 300;
	int poissonDepth = 11;
	const int alveoli = 8 * 30000;
	simulation *s = new simulation();
	status * st = new status();
	try {
		s->extendBronchialTreeV2(
			buildVolumes,
			buildCenterline,
			true,
			build1DModel,
			surfaceSampling,
			buildNormals,
			build3DModel,
			refinements,
			segmentation,
			depth,
			volumeDepth,
			density,
			poissonDepth,
			path,
			boundaryR,
			boundaryL,
			trachea,
			existingModel,
			extractedCenterlineModel,
			st,
			true,
			leftLungModel,
			rightLungModel
		);
	}
	catch (const std::exception& e) {
		std::cout << " a standard exception was caught, with message '" << e.what() << "'\n";
	}
	delete s;
```











