#include"Process.h"
#include<string>
#include<iostream>
#include<fstream>
#include<io.h>
using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		cerr << "Syntax: PP_Code.exe source_mesh.ovm initial_mesh.ovm" << endl;
		return -1;
	}
	if (3 == argc)
	{
		const string source_mesh = argv[1];
		const string initial_mesh = argv[2];

		Process pro(source_mesh, initial_mesh);

		cout << " foldover-free volumetric mapping construction begin ..." << endl;
		pro.pck();
		cout << " foldover-free volumetric mapping construction finish !" << endl;

		//system("pause");
		return 0;
	}
}
