#include <stdio.h>
#include <vector>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "Vec3.h"

#ifndef FILEREADER
#define FILEREADER


class fileReader
{
    static void checkFileExt(const std::string inFile, const std::string ext)
    {
        //convert file name to string.
        //old code took a char array as an input
        //so the string stream is still kept here;
        std::stringstream ss;
        std::string fileName;
        ss << inFile;
        ss >> fileName;

        //get file extension
        int period = fileName.find_last_of(".");
        if(period < 0)
            throw "Error: No file extension given\n";

        std::string fileExt = fileName.substr(period + 1);
        std::transform(fileExt.begin(), fileExt.end(), fileExt.begin(), ::tolower);

        if (fileExt.compare(ext.c_str()))
            throw "Error: Incorrect file type, must be XYZ\n";
    };

    public:
    typedef std::vector<std::vector<std::string> > flines_p;
    typedef std::vector<std::string > flines;
    static std::vector<std::vector<std::string> > readFile(const std::string fileName)
    {
        
        //using namespace std;
        using std::ifstream;
        using std::vector;
        using std::string;

        ifstream inFile;
        string inBuffer, token;

        vector<string> lineContents;
        vector<vector<string> > saveContents;
        std::stringstream streamBuffer;

        // read in file line by line
        inFile.open(fileName.c_str());
        if((bool)inFile)
        {   
            while(getline(inFile, inBuffer))
            {
                // convert string to buffer and parse by spaces
                streamBuffer << inBuffer;
                while(streamBuffer >> token)
                {
                    
                    lineContents.push_back(token);
                }
                saveContents.push_back(lineContents);

                // clear buffers and streams
                lineContents.clear();
                streamBuffer.clear();
            }

            inFile.close();
            return saveContents;
        }
        else
            throw string("Could not open file");
    }

    static std::vector<std::string > readFile_unparsed(const std::string fileName)
    {
        //using namespace std;

        std::ifstream inFile;
        std::string inBuffer;

        std::vector<std::string > saveContents;

        // read in file line by line
        inFile.open(fileName.c_str());
        while(getline(inFile, inBuffer))
        {
            saveContents.push_back(inBuffer);
        }

        inFile.close();
        return saveContents;
    }

    static void readXYZ(const std::string fileName, std::vector<Vec3> &coords, std::vector<std::string> &atoms)
    {
        std::ifstream molFile;
        std::vector<double> tmpVec;
        std::string type;
        double x, y, z;
        int atomIdx = 0;
        coords.clear();
        atoms.clear();

        try
        {
        	checkFileExt(fileName, "xyz");
        }
        catch (const char* msg)
        {
        	std::cerr << msg;
        	throw "Error: Could not import molecule file\n";
        }
        //open molecule file for reading
        molFile.open(fileName.c_str());

        // skip first two lines
        std::string tmpStr;
        getline(molFile, tmpStr);
        getline(molFile, tmpStr);

        //read in molecule file
        while (molFile >> type >> x >> y >> z)
        {
            Vec3 pos(x, y, z);
            coords.push_back(pos);
            atoms.push_back(type);
        }
        molFile.close();
    }
};

#endif