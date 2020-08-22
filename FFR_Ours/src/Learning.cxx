#include "Learning.h"
#include "LearningImpl.h"

Learning::Learning()
{
	limpl = new LearningImpl;
}
Learning::~Learning()
{
	delete limpl;
}
