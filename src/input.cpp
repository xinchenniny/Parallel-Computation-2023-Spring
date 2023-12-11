#include  "input.h"
#include  "timer.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
Input::Input(string filename) {
    timer::tick("INPUT","Input");
    parseFile(filename);
    timer::tick("INPUT","Input");
}

void Input::parseFile(string filename) {
    timer::tick("INPUT","parseFile");
    ifstream file(filename);

    // 如果文件不能打开，输出错误信息
    if(!file.is_open()) {
        cerr << "Error: Cannot open input file!" << endl;
        exit(1);
    }

    map<string, string> data; // 存储键值对

    string line;
    while(getline(file, line)) {
        if(line.empty() || line[0] == '#') // 忽略空白行和注释行
            continue;
        size_t pos = line.find(' ');
        string key = line.substr(0, pos); // 获取键名
        string value = line.substr(pos+1); // 获取键值
        data[key] = value; // 添加到键值对中
    }

    // 从键值对中读取各个参数的数值
    isHexahedral = stoi(data["isHexahedral"]);
    lx = stod(data["lx"]);
    ly = stod(data["ly"]);
    lz = stod(data["lz"]);
    thetaxy = stod(data["thetaxy"]);
    thetayz = stod(data["thetayz"]);
    thetaxz = stod(data["thetaxz"]);
    support_SH = stoi(data["support_SH"]);
    diago_lib = data["diago_lib"];
    support_Periodic_Boundary = stoi(data["support_Periodic_Boundary"]);
    multi_parallel_strategies = stoi(data["multi_parallel_strategies"]);
    points_path = data["points_path"];
    venergy_path = data["v_path"];
    distribution_path = data["distribution_path"];

    file.close();
    timer::tick("INPUT","parseFile");

}

std::vector<Input::Point> Input::readpointsFile(const std::string& filename, int limit) {
    timer::tick("INPUT","readpointsFile");
    std::vector<Point> points;
    std::ifstream file(filename);
    std::string line;
    points.clear();  // 初始化为空
    points.reserve(limit);  // 预分配limit个元素的空间    
    while (std::getline(file, line)) {
        if (limit > 0 && count >= limit) break;  // 如果给定了上限，则最多读取limit个点
        std::stringstream ss(line);
        double x, y, z;
        char comma;
        Point point; 
        ss >> comma >> point.x >> comma >> point.y >> comma >> point.z >> comma;
        points.push_back(point);      
        std::cout << "(" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
        count++;
    }
    cout << "Read in " << count <<" points."<<endl;
    timer::tick("INPUT","readpointsFile");
    return points;
}

void  Input::readDirstibution(const std::string& filename){
    timer::tick("INPUT","readDirstibution");
   std::ifstream file(filename, std::ios::in);
        if (!file) std::cerr << "Can't open the file\n";
        std::string line;

        // 读取文本中的数据并存入 vector<double> 中
        while (getline(file, line)) {
            if (line.empty() || line[0] == '#') continue; // 剔除空行和注释行

            std::stringstream sstream(line);
            if (line.find("cutoff") != string::npos) {
                sscanf(line.c_str(), "cutoff %lf", &cutoff);
            } else if (line.find("dr") != string::npos) {
                sscanf(line.c_str(), "dr %lf", &dr);
            } else if (line.find("mesh") != string::npos) {
                sscanf(line.c_str(), "mesh %d", &mesh);
            } else if (line.find("l") != string::npos) {
                sscanf(line.c_str(), "l %d", &l);
            } else {
                double temp;
                while (!sstream.eof()) {
                    sstream >> temp;
                    if (sstream.fail()) break; // 处理数据读取失败的情况
                    f.push_back(temp); // 存储径向分布函数至 vector<double> f 中
                    char ch = sstream.peek(); // 查看下一个字符是否为分隔符（逗号或空格）
                    if (ch == ',' || ch == ' ') sstream.ignore();
                }
            }
        }
        file.close();
       timer::tick("INPUT","readDirstibution");
}

void Input::readVfile(const std::string& filename){
       timer::tick("INPUT","readVfile");
    ifstream infile(filename);
    string str;

    while (infile >> str) {
        if (str == "nx") {
            infile >> nx;
        }
        else if (str == "ny") {
            infile >> ny;
        }
        else if (str == "nz") {
            infile >> nz;
        }
        else if (str == "V:") {
            for (int i = 0; i < nx * ny * nz; i++) {
                double value;
                infile >> value;
                V.push_back(value);
            }
        }
    }

    // 输出读取得到的 V 数组
    for (int i = 0; i < nx * ny * nz; i++) {
        //cout << V[i] << " ";
    }
    cout << endl;
    timer::tick("INPUT","readVfile");

}