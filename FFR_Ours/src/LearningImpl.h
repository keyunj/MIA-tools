#ifndef __LearningImpl_h__
#define __LearningImpl_h__

#include <vector>
#include <string>
#include "opencv2/core/core.hpp"
#include "opencv2/ml/ml.hpp"

class LearningImpl
{
public:
	LearningImpl();
	~LearningImpl();
	bool LoadLandmarkClassifiers(int num, std::string& workPath);
	bool LoadLumenWallClassifiers(std::string& workPath);
	bool LoadBoundaryClassifiers(int num);
	std::vector<CvBoost> lmBoost;
	CvBoost				 lwBoost;
	std::vector<CvBoost> bdBoost;
private:
	bool lmLoaded;
	bool lwLoaded;
	bool bdLoaded;
};


#endif