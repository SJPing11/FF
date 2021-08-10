#include"Parafun.h"
#include<string>
#include<iostream>
#include<fstream>
#include<io.h>
using namespace std;

int main(int argc, char *argv[])
{
#if 1
	if (argc < 3)
	{
		cerr << "Syntax: " << argv[0] << " <input mesh> <initial mesh> <dim = 2>" << endl;
		return -1;
	}
	if (3 == argc)
	{
		const string input_mesh = argv[1];
		const string initial_mesh = argv[2];

		Parafun pro(input_mesh, initial_mesh, 2);

		cout << input_mesh << " plane process projective begin ..." << endl;
		pro.pck();
		cout << input_mesh << " plane process projective finish !" << endl;

		system("pause");
		return 0;
	}
	const string input_mesh = argv[1];
	const string initial_mesh = argv[2];

	Parafun pro(input_mesh, initial_mesh, int(argv[3]));

	cout << input_mesh << " surface process projective begin ......" << endl;
	pro.pck();
	cout << input_mesh << " surface process projective finish !" << endl;

	system("pause");
	return 0;
#else
	//string path = "E:/Disk_Topology_Mesh/Mesh1000/";
	string path = "E:/Code/AKVFParam/";
	string pathsearch = path + "*.obj";
	string save_dir_str = "E:/ceshi/";
	vector<string> filesname;

	_finddata_t fileinfo;
	auto file_h = _findfirst(pathsearch.c_str(), &fileinfo);
	do
	{
		filesname.push_back(path + fileinfo.name);
	} while (_findnext(file_h, &fileinfo) == 0);

	cout << "there exists " << filesname.size() << " models need to be processed;" << endl;

	for (size_t i_file = 0; i_file < filesname.size(); i_file++)
	{
		//string modelname = filesname[i_file].substr(17, 5);
		cout << filesname[i_file] << "---------------------------------begin......." << endl;
		Parafun parafun_m(filesname[i_file], save_dir_str);
		parafun_m.run_composite();

		//parafun_m.run_withoutintp();
		//parafun_m.run_init_dis_distribut();
		//parafun_m.run_adaptive();
		cout << filesname[i_file] << "---------------------------------finish......." << endl;
	}
	system("pause");
#endif // 0


}
