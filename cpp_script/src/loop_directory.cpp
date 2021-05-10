#include <iostream>
#include <vector>
#include <string>
#include <experimental/filesystem>

using namespace std;

int main()
{
    string path = "../../python_script/WT_R164S_65_213/pdb/WT_offset_0/";
    for (const auto & file : experimental::filesystem::directory_iterator(path))
        cout << file.path() << endl;
    return EXIT_SUCCESS;
}



