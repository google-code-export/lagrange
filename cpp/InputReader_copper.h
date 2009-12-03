/*
 * InputReader.h
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 */

#ifndef INPUTREADER_COPPER_H_
#define INPUTREADER_COPPER_H_

#include <vector>
#include <string>
using namespace std;


class InputReader{
	public:
	InputReader();
		vector<TreeTemplate<Node> *> readMultipleTreeFile(string filename);
		map<string,vector<int> > readStandardInputData(string filename);
		void checkData(map<string,vector<int> >,vector<TreeTemplate<Node> *>);
		int nareas;
		int nspecies;
};

#endif /* INPUTREADER_H_ */
