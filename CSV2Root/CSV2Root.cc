#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TString.h>

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t\n\r")
{
	const auto strBegin = str.find_first_not_of(whitespace);
	if (strBegin == std::string::npos) return ""; // no content

	const auto strEnd = str.find_last_not_of(whitespace);
	const auto strRange = strEnd - strBegin + 1;

	return str.substr(strBegin, strRange);
}

std::string reduce(const std::string& str,
                   const std::string& fill = " ",
                   const std::string& whitespace = " \t\n\r")
{
	// trim first
	auto result = trim(str, whitespace);

	// replace sub ranges
	auto beginSpace = result.find_first_of(whitespace);
	while (beginSpace != std::string::npos) {
		const auto endSpace = result.find_first_not_of(whitespace, beginSpace);
		const auto range = endSpace - beginSpace;

		result.replace(beginSpace, range, fill);

		const auto newStart = beginSpace + fill.length();
		beginSpace = result.find_first_of(whitespace, newStart);
	}

	return result;
}

std::vector<std::string> tokenize(const std::string& line,
                                  const char delim)
{
	std::vector<std::string> tokens;

	std::istringstream ss(line);
	std::string item;
	while (std::getline(ss, item, delim)) {
		item = reduce(item, "_");
		tokens.push_back(item);
	}
	return tokens;
}

int main(int argc, char* argv[])
{
	if (argc!=4) {
		std::cerr << "usage: CSV2Root <line number with column names> <first line with data> <filename>" << std::endl;
		return -1;
	}
	
	int nameline = std::atoi(argv[1]);
	int firstdataline = std::atoi(argv[2]);
	std::ifstream ifile(argv[3]);
	
	std::string line;
	int linenum = 0;
	int colnum = 0;
	
	TFile * ofile = new TFile("output.root", "RECREATE");
	TTree * otree = new TTree("CSV2Root", "CSV2Root");
	
	std::map<int,std::string> branchnames;
	double* coldata; 
	std::map<int,TBranch*> branches;
	
	while (std::getline(ifile, line)) {
		linenum++;
		
		if (linenum==nameline) {
			std::vector<std::string> tokens = tokenize(line, ';');
			colnum = 0;
			for (std::vector<std::string>::iterator it=tokens.begin();
			     it!=tokens.end();
			    ++it) {
					branchnames[colnum] = *it;
					colnum++;
			}
			coldata = new double[colnum];
			colnum = 0;
			for (std::vector<std::string>::iterator it=tokens.begin();
			     it!=tokens.end();
			    ++it) {
			    TBranch * branch = otree->Branch((*it).c_str(), &coldata[colnum], Form("%s/D", (*it).c_str()));
					branches[colnum] = branch;
					colnum++;
			}
		}
		
		if (linenum>=firstdataline) {
			std::vector<std::string> tokens = tokenize(line, ';');
			colnum = 0;
			for (std::vector<std::string>::iterator it=tokens.begin();
			     it!=tokens.end();
			    ++it) {
			  coldata[colnum] = std::atof((*it).c_str());
			  colnum++;
			}
			otree->Fill();
		}
	}
	
	ofile->Write();
	ofile->Close();
	delete ofile;
	
	delete [] coldata;

	return 1;
}