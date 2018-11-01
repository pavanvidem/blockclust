/*
 * Read.h
 *
 *  Created on: Nov 9, 2012
 *      Author: videmp
 */

#ifndef TAG_H_
#define TAG_H_
#include <string>
using namespace std;

class Read{
	public:
		Read();
		Read(string rdId, unsigned long rdStart, unsigned long rdEnd, double rdExpression, double noOfMappings, double readCount);
		virtual ~Read();
		string getId();
		unsigned long getStart();
		unsigned long getEnd();
		unsigned short getLength();
		double getExpression();
		double getNumberOfMappings();
		double getReadCount();
		double getArea();
	private:
		string readId;
		unsigned long readStart;
		unsigned long readEnd;
		double readExpression;
		double numberOfMappings;
		double readCount;
};

#endif /* TAG_H_ */
