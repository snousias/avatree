#include "lungModelling.h"
#include <iostream>

#include "itkImage.h"
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkMetaImageIO.h>
#include <itkImageRegionIterator.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/selection.h>
#include <fstream>



int main(void) {
	std::string root = "F:\\_Datasets\\CT-MRI-DATABASE\\_VESSEL\\Scans\\VESSEL12_08\\";
	std::string path = root;
	std::string mhdLungFileName = root + "VESSEL12_08.mhd";
	std::string rawLungFileName = root + "VESSEL12_08.raw";
	std::string mhdLungMaskFileName = root + "VESSEL12_08_Lungs.mhd";
	std::string rawLungMaskFileName = root + "VESSEL12_08_Lungs.raw";
	int seed, label, lim;
	std::vector<int> labelsArray = { 100,200 };

	lim = 500;

	bool buildVolumes = true;
	bool buildCenterline = true;
	bool build1DModel = true;
	bool surfaceSampling = true;
	bool buildNormals = true;
	bool build3DModel = true;
	bool refinements = false;
	bool segmentation = false;
	int depth = 12;
	int volumeDepth = 16;
	int density = 300;
	int poissonDepth = 11;
	const int alveoli = 8 * 30000;

	typedef itk::ImageIOBase::IOComponentType ScalarPixelType;
	voxelSpace * Vox = new voxelSpace();
	itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
		mhdLungMaskFileName.c_str(), itk::ImageIOFactory::ReadMode);
	if (!imageIO)
	{
		std::cerr << "Could not CreateImageIO for: " << mhdLungMaskFileName << std::endl;
		return EXIT_FAILURE;
	}
	imageIO->SetFileName(mhdLungMaskFileName.c_str());
	imageIO->ReadImageInformation();
	std::cout << "Pixel Type is " << imageIO->GetComponentType() << std::endl;
	const size_t numDimensions = imageIO->GetNumberOfDimensions();
	std::cout << "Number of Dimensions: " << numDimensions << std::endl;
	std::cout << "Component size: " << imageIO->GetComponentSize() << std::endl;

	std::cout << "Pixel type: " << imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) << std::endl;

	std::cout << "Byte Order: " << imageIO->GetByteOrder() << std::endl;
	std::cout << "Dimension X " << imageIO->GetDimensions(0) << std::endl;
	std::cout << "Dimension Y " << imageIO->GetDimensions(1) << std::endl;
	std::cout << "Dimension Z " << imageIO->GetDimensions(2) << std::endl;

	std::cout << "Dimension X " << imageIO->GetSpacing(0) << std::endl;
	std::cout << "Dimension Y " << imageIO->GetSpacing(1) << std::endl;
	std::cout << "Dimension Z " << imageIO->GetSpacing(2) << std::endl;

	std::cout << "Component type: " << imageIO->GetComponentTypeAsString(imageIO->GetComponentType()) << std::endl;

	std::cout << "Byte Order: " << imageIO->GetByteOrderAsString(imageIO->GetByteOrder()) << std::endl;

	float dx = imageIO->GetSpacing(0);
	float dy = imageIO->GetSpacing(1);
	float dz = imageIO->GetSpacing(2);
	int LX = imageIO->GetDimensions(0);
	int LY = imageIO->GetDimensions(1);
	int LZ = imageIO->GetDimensions(2);

	std::vector<std::vector<std::vector<int>>> A;
	Vox->resize3DMat(A, mDiscretization(LZ, LY, LX));

	VectorXf  a2 = VectorXf::Zero(imageIO->GetDimensions(0)*imageIO->GetDimensions(1)*imageIO->GetDimensions(2));
	std::ifstream myData(rawLungMaskFileName, std::ios::binary);
	unsigned char value;
	int i = 0;
	char buf[sizeof(unsigned char)];
	while (myData.read(buf, sizeof(buf)))
	{
		memcpy(&value, buf, sizeof(value));
		//std::cout << value << " ";
		a2(i) = value;
		i++;
	}
	std::cout << std::endl << "Total count: " << i << std::endl;
	int idx = 0;
	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			for (int k = 0; k < A[i][j].size(); k++) {
				int val = a2(idx++);
				A[i][j][k] = val;
			}
		}
	}

	std::vector<mCellSimplified> Asparse;
	Vox->Matrix2SparseRepresentationByValue(A, Asparse, 1.0);

	label = labelsArray[0];
	std::vector<mIndex> toInvestigate;
	seed = Vox->findSeed(A, Asparse, 1);
	mIndex seedIndex = Asparse[seed].index;
	A[seedIndex.x][seedIndex.y][seedIndex.z] = label;

	toInvestigate.clear();
	toInvestigate.push_back(seedIndex);
	int regionGrowingFront = 0;
	for (int pp = 0; pp < toInvestigate.size(); pp++) {
		mIndex * ind = &toInvestigate[pp];
		if (ind->isInsideMatrix(A, 2)) {
			bool found = Vox->regionHasNotLabel(ind, A, 2, 0);
			bool isNeighb = Vox->regionHasOnlyLabel(ind, A, 1, label);
			if (found && isNeighb) {
				A[ind->x][ind->y][ind->z] = label;
			}
		}
		if (pp == toInvestigate.size() - 1) {
			regionGrowingFront++;
			int c = regionGrowingFront;
			int r1, r2, r3;
			toInvestigate.clear();
			pp = 0;
			for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1 = r1 + 2 * c) {
				for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2++) {
					for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
			for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2 = r2 + 2 * c) {
				for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1++) {
					for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
			for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3 = r3 + 2 * c) {
				for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1++) {
					for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
		}
		if (regionGrowingFront == lim) {
			regionGrowingFront = 0;
			toInvestigate.clear();
		}
	}

	label = labelsArray[1];

	seed = Vox->findSeed(A, Asparse, 1);
	seedIndex = Asparse[seed].index;
	A[seedIndex.x][seedIndex.y][seedIndex.z] = label;
	toInvestigate.clear();
	toInvestigate.push_back(seedIndex);
	regionGrowingFront = 0;

	for (int pp = 0; pp < toInvestigate.size(); pp++) {
		mIndex * ind = &toInvestigate[pp];
		if (ind->isInsideMatrix(A, 2)) {
			if (Vox->regionHasNotLabel(ind, A, 2, 0) && Vox->regionHasNotLabel(ind, A, 2, label) && Vox->regionHasOnlyLabel(ind, A, 1, label)) {
				A[ind->x][ind->y][ind->z] = label;
			}
		}
		if (pp == toInvestigate.size() - 1) {
			regionGrowingFront++;
			int c = regionGrowingFront;
			int r1, r2, r3;
			toInvestigate.clear();
			pp = 0;
			for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1 = r1 + 2 * c) {
				for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2++) {
					for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
			for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2 = r2 + 2 * c) {
				for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1++) {
					for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
			for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3 = r3 + 2 * c) {
				for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1++) {
					for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
		}
		if (regionGrowingFront == lim) {
			regionGrowingFront = 0;
			toInvestigate.clear();
		}
	}

	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			for (int k = 0; k < A[i][j].size(); k++) {
				if (A[i][j][k] == 1) {
					A[i][j][k] = 0;
				}
			}
		}
	}

	std::cout << A.size() << " " << A[0].size() << " " << A[0][0].size() << std::endl;

	volume * result2 = new volume();
	while (result2->vertices.size() < alveoli) {
		float XSampling = result2->functions.uniformRandom()  * (LX - 1);
		float YSampling = result2->functions.uniformRandom()  * (LY - 1);
		float ZSampling = result2->functions.uniformRandom()  * (LZ - 1);
		//std::cout << XSampling << " " << YSampling << " " << ZSampling << std::endl;
		if (A[floor(ZSampling)][floor(YSampling)][floor(XSampling)] == labelsArray[0]) {
			result2->vertices.push_back({ XSampling*dx,YSampling*dy,ZSampling*dz });
		}
	}
	/*for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			for (int k = 0; k < A[i][j].size(); k++) {
				if (A[i][j][k] == labelsArray[0]) {
					result2->vertices.push_back({ (float)k*dx,(float)j*dy,(float)i*dz });
				}
			}
		}
	}*/
	result2->exportToFile("LL.obj");

	volume * result3 = new volume();
	while (result3->vertices.size() < alveoli) {
		float XSampling = result3->functions.uniformRandom() * (LX - 1);
		float YSampling = result3->functions.uniformRandom() * (LY - 1);
		float ZSampling = result3->functions.uniformRandom() * (LZ - 1);
		//std::cout << XSampling <<  " " <<YSampling << " "  << ZSampling << std::endl;
		if (A[floor(ZSampling)][floor(YSampling)][floor(XSampling)] == labelsArray[1]) {
			result3->vertices.push_back({ XSampling*dx,YSampling*dy, ZSampling*dz });
		}
	}
	/*for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			for (int k = 0; k < A[i][j].size(); k++) {
				if (A[i][j][k] == labelsArray[1]) {
					result3->vertices.push_back({ (float)k*dx,(float)j*dy,(float)i*dz });
				}
			}
		}
	}*/
	result3->exportToFile("RL.obj");

	//std::string boundaryMeshL = path + "VESSEL12_08_Lungs_L_3.obj";
	//std::string boundaryMeshR = path + "VESSEL12_08_Lungs_R_3.obj";

	std::string existingModelMesh = path + "existingModel2.obj";
	std::string extractedCenterline = path + "VESSEL12_08_Centerline.vtk";

	std::cout << std::endl << "Model :: Right host mesh " << std::endl;
	dotObj * hostR = new dotObj();
	std::cout << std::endl << "Model :: Right boundary " << std::endl;
	dotObj * boundaryR = new dotObj();
	//boundaryR->initializeFromFile(boundaryMeshR);
	std::cout << std::endl << "Model :: Left host mesh " << std::endl;
	dotObj * hostL = new dotObj();
	std::cout << std::endl << "Model :: Left boundary " << std::endl;
	dotObj * boundaryL = new dotObj();
	//boundaryL->initializeFromFile(boundaryMeshL);
	std::cout << std::endl << "Model :: Existing geometry " << std::endl;
	dotObj *existingModel = new dotObj();
	existingModel->initializeFromFile(existingModelMesh);

	dotObj * extractedCenterlineModel = new dotObj();
	//extractedCenterlineModel->initializeFromVTK(extractedCenterline);
	*extractedCenterlineModel = existingModel->mcfskel();
	extractedCenterlineModel->exportToFile(path + "cccline");

	dotObj * trachea = new dotObj();
	if (strlen(path.c_str()) > 2) {
	}
	else {
		std::cout << "Error: Enter valid workspace path" << std::endl;
		return 0;
	}
	simulation *s = new simulation();
	status * st = new status();
	try {
		//s->extendbronchialtreev2(
		//	buildvolumes,
		//	buildcenterline,
		//	build1dmodel,
		//	surfacesampling,
		//	buildnormals,
		//	build3dmodel,
		//	refinements,
		//	segmentation,
		//	depth,
		//	volumedepth,
		//	density,
		//	poissondepth,
		//	path,
		//	boundaryr,
		//	boundaryl,
		//	trachea,
		//	existingmodel,
		//	st,
		//	true);

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
			result2,
			result3
		);
	}
	catch (const std::exception& e) {
		std::cout << " a standard exception was caught, with message '" << e.what() << "'\n";
	}
	delete s;

	return 0;

	dotObj * fullTreeExtendedModel = new dotObj();
	dotObj * rec = new dotObj();

	//fullTreeExtendedModel->initializeFromFile(path + "_7_fullTreeOversampled.obj");
	fullTreeExtendedModel->initializeFromFile(path + "_7_fullTree.obj");
	//fullTreeExtendedModel = extractedCenterlineModel;
	fullTreeExtendedModel->scale(1 / dz, 1 / dy, 1 / dx);
	fullTreeExtendedModel->graphBased1DModelAnalysis();
	//fullTreeExtendedModel->mAnalysis.exportGraphFeatures(path+"distributionOfAnglesnLengths.csv");

	fullTreeExtendedModel->mAnalysis.graph.ggraph2Model(rec);
	rec->exportToFile(path + "_7_rec");
	std::cout << "Grid Generation complete" << std::endl;
	for (int i = 0; i < fullTreeExtendedModel->mAnalysis.graph.nodes.size(); i++) {
		Vector3f * p = &fullTreeExtendedModel->mAnalysis.graph.nodes[i].position;
		int iz = int(p->x());
		int iy = int(p->y());
		int ix = int(p->z());
		if ((A[ix][iy][iz] == labelsArray[0]) || (A[ix][iy][iz] == labelsArray[1])) {
			A[ix][iy][iz] = fullTreeExtendedModel->mAnalysis.graph.nodes[i].generation;
		}
	}
	std::cout << "Distribution complete" << std::endl;

	std::vector<mCellSimplified> Bsparse;
	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			for (int k = 0; k < A[i][j].size(); k++) {
				if ((A[i][j][k] > 0) && (A[i][j][k] < 20)) {
					Bsparse.push_back(mCellSimplified(i, j, k, A[i][j][k]));
				}
			}
		}
	}
	std::vector<std::vector<std::vector<pCell>>> Dist2Gen;
	Vox->resize3DMat(Dist2Gen, mDiscretization(LZ, LY, LX));

	int gstart = 4;
	int gend = 9;
	for (int gensearch = gstart; gensearch < gend; gensearch++)
	{
		for (int i = 0; i < A.size(); i++) {
			for (int j = 0; j < A[i].size(); j++) {
				for (int k = 0; k < A[i][j].size(); k++) {
					if ((A[i][j][k] > 0)) {
						float minDist = 1000;
						for (int g = 0; g < Bsparse.size(); g++) {
							if (Bsparse[g].value == gensearch) {
								Vector3f m1 = Vector3f(float(i), float(j), float(k));
								Vector3f m2 = Vector3f(float(Bsparse[g].index.x), float(Bsparse[g].index.y), float(Bsparse[g].index.z));
								float cDist = (m1 - m2).norm();
								if (cDist < minDist) {
									minDist = cDist;
								}
							}
						}
						//std::cout << minDist << std::endl;
						//Dist2Gen[i][j][k] = minDist;
						float distributionAmplitude = 200.0;
						float pdf_gaussian = (1 / (10 * sqrt(2 * M_PI))) * exp(-0.5 * pow((minDist - 0) / 10, 2.0));
						Dist2Gen[i][j][k].properties[gensearch].value = distributionAmplitude * pdf_gaussian;
						//Dist2Gen[i][j][k].properties[1].value = distributionAmplitude;
					}
					else {
						Dist2Gen[i][j][k].properties[gensearch].value = 0;
						//Dist2Gen[i][j][k].properties[1].value = 0;
					}
				}
			}
		}
	}

	for (unsigned int g = gstart; g < gend; g++) {
		short * buffer = new short[LX*LY*LZ];
		int ccc = 0;
		for (int i = 0; i < Dist2Gen.size(); i++) {
			for (int j = 0; j < Dist2Gen[i].size(); j++) {
				for (int k = 0; k < Dist2Gen[i][j].size(); k++) {
					buffer[ccc] = (short)Dist2Gen[i][j][k].properties[g].value;
					ccc++;
				}
			}
		}
		std::string fnmfinal = path + "probability_map_" + std::to_string(g) + ".raw";
		fullTreeExtendedModel->functions.write_file_binary(fnmfinal, reinterpret_cast<char const *>(buffer), sizeof(short) * LX*LY*LZ);
		delete[] buffer;
	}

	return 0;
}

