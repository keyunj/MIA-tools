#include "LearningImpl.h"

#include <sstream>
#include <Windows.h>

LearningImpl::LearningImpl()
{
	lmLoaded = false;
	lwLoaded = false;
	bdLoaded = false;
}

LearningImpl::~LearningImpl()
{
}

bool LearningImpl::LoadLandmarkClassifiers(int num, std::string& workPath)
{
	if(lmLoaded) return true;

	std::ostringstream strstm; 
	lmBoost.resize(num);

	cv::FileStorage fs(workPath + "\\lvcorlmclassifier.yml.gz", cv::FileStorage::READ);
	bool temp = fs.isOpened();
	for(int id=0; id<num; id++)
	{
		strstm.str("");
		strstm << "Classifier_" << id;
		lmBoost[id].read(*fs, *fs[strstm.str().c_str()]);
		if( !lmBoost[id].get_weak_predictors() )
		{
			std::cerr << "Could not read " << strstm.str() << std::endl;
			return false;
		}
	}
	fs.release();
	lmLoaded = true;
	
	return true;
}

bool LearningImpl::LoadLumenWallClassifiers(std::string& workPath)
{
	if(lwLoaded) return true;

	cv::FileStorage fs(workPath + "\\lumenwallclassifier.yml.gz", cv::FileStorage::READ);
	{
		std::ostringstream strstm; 
		strstm.str("");
		strstm << "Classifier";
		lwBoost.clear();
		lwBoost.read(*fs, *fs[strstm.str().c_str()]);
		if( !lwBoost.get_weak_predictors() )
		{
			std::cerr << "Could not read " << strstm.str() << std::endl;
			return false;
		}
	}
	fs.release();
	
	lwLoaded = true;
	return true;
}

bool LearningImpl::LoadBoundaryClassifiers(int num)
{
	if(bdLoaded) return true;

	std::ostringstream strstm; 
	bdBoost.resize(num);

	cv::FileStorage fs("lvbdclassifier.yml.gz", cv::FileStorage::READ);
	for(int id=0; id<num; id++)
	{
		strstm.str("");
		strstm << "Classifier_" << id;
		bdBoost[id].read(*fs, *fs[strstm.str().c_str()]);
		if( !bdBoost[id].get_weak_predictors() )
		{
			std::cerr << "Could not read " << strstm.str() << std::endl;
			return false;
		}
	}
	fs.release();

	bdLoaded = true;
	return true;
}
