// TestingFunc.cpp: implementation of the TestingFunc class.
//
//////////////////////////////////////////////////////////////////////
//#include <iostream>
//#include <fstream>
//#include <string>

//using namespace std;
#include "StrVis.hpp"

string convertToZeroPaddedString(int number, int width){
    stringstream ss;
    ss << setw(width) << setfill('0') << number;
    return ss.str();
}

void StrVis(string exepath, string pdbFilePath1, string pdbFilePath2, string pdboutFilePath, int imageid){
//int main() {
    // PDB形式の座標データを含むファイルのパス
    //std::string pdbFilePath1 = "/Users/tsukiyamashou/Documents/Project_RNA3D/RNARefiner/globalRNArefiner/src/ex/out/ss_tmp.pdb";
    //std::string pdbFilePath2 = "/Users/tsukiyamashou/Documents/Project_RNA3D/RNARefiner/globalRNArefiner/src/ex/out/ss_init.pdb";
    //std::string pdboutFilePath =  "/Users/tsukiyamashou/Documents/Project_RNA3D/RNARefiner/globalRNArefiner/vissrc/current_str.png";
    // PyMOLに接続するためのスクリプトを生成
    std::ofstream pymolScript(exepath+"/pymol_script.py");
    pymolScript << "from pymol import cmd\n";
    pymolScript << "import warnings\n";
    pymolScript << "warnings.simplefilter('ignore')\n";
    pymolScript << "cmd.load('" << pdbFilePath1 << "')\n";
    pymolScript << "cmd.load('" << pdbFilePath2 << "')\n";
    pymolScript << "cmd.show('sticks', 'decoy')\n";
    pymolScript << "cmd.show('sticks', 'init_decoy')\n";
    pymolScript << "cmd.align('decoy', 'init_decoy')\n";
    pymolScript << "cmd.zoom()\n";
    //pymolScript << "cmd.show('sticks')\n";
    //pymolScript << "cmd.color('green', 'polymer')\n";
    pymolScript << "cmd.color('red', 'decoy')\n";
    pymolScript << "cmd.color('blue', 'init_decoy')\n";
    pymolScript << "cmd.hide('cartoon', 'all')\n";
    
    //pymolScript << "cmd.show('surface')\n";
    pymolScript << "cmd.bg_color('white')\n";
    pymolScript << "cmd.viewport(400, 400)\n";
    pymolScript << "cmd.refresh()\n";

    pymolScript << "cmd.png('" << pdboutFilePath + "/refined_str_" + convertToZeroPaddedString(imageid, 7) << ".png', width=400, height=400, dpi=300)\n";  // 画像保存コマンドを追加
    pymolScript.close();

    // PyMOLを起動してスクリプトを実行
    std::string command = "pymol -cq " + exepath + "/pymol_script.py";
    int result = system(command.c_str());

    if (result != 0) {
        std::cerr << "Failed to execute PyMOL script." << std::endl;
        exit(EXIT_FAILURE);  // エラーで終了する
    }

    //return 0;
}



/*
void StrVis(string exepath){
    std::string command = "pymol -cq " + exepath + "/pymol_script.py";
    int result = system(command.c_str());

    if (result != 0) {
        std::cerr << "Failed to execute PyMOL script." << std::endl;
        exit(EXIT_FAILURE);  // エラーで終了する
    }
}
*/
	
/*
int main(){
    std::string pdbFilePath1 = "/Users/tsukiyamashou/Documents/Project_RNA3D/RNARefiner/globalRNArefiner/src/ex/out/ss_tmp.pdb";
    std::string pdbFilePath2 = "/Users/tsukiyamashou/Documents/Project_RNA3D/RNARefiner/globalRNArefiner/src/ex/out/ss_init.pdb";
    std::string pdboutFilePath =  "/Users/tsukiyamashou/Documents/Project_RNA3D/RNARefiner/globalRNArefiner/vissrc";
    int imageid=0;

    StrVis(pdbFilePath1, pdbFilePath2, pdboutFilePath, imageid);
    return 0;
}
*/

