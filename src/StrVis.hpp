// TestingFunc.h: interface for the TestingFunc class.
//
//////////////////////////////////////////////////////////////////////

//#if !defined(AFX_TESTINGFUNC_H__5ADC12DB_C425_45A5_A3B8_7DBBC34CAE4A__INCLUDED_)
//#define AFX_TESTINGFUNC_H__5ADC12DB_C425_45A5_A3B8_7DBBC34CAE4A__INCLUDED_
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> // setw
#include <sstream> // stringstream
using namespace std;

#ifndef STRVIS
#define STRVIS

string convertToZeroPaddedString(int number, int width);
void StrVis(string exepath, string pdbFilePath1, string pdbFilePath2, string pdboutFilePath, int imageid);
//void StrVis(string exepath);

#endif

