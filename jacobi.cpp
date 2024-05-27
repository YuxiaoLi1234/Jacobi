#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>

std::set<std::pair<int, int> > false_edges;
std::vector<std::vector<int> > adjacency;
std::map<std::pair<int, int>, int> is_critical;
std::map<int, int> frequency_map;

std::map<std::pair<int, int>, int> decp_is_critical;
std::vector<double> input_f;
std::vector<double> input_g;
std::vector<double> decp_f;
std::vector<double> decp_g;
std::vector<double> getdata(std::string filename ){
    
    std::ifstream inputFile(filename, std::ios::binary);

    // 检查文件是否成功打开
    if (!inputFile) {
        std::cerr << "无法打开文件 example.bin" << std::endl;

    }

    // 移动文件指针到文件末尾以确定文件大小
    inputFile.seekg(0, std::ios::end);
    std::streampos fileSize = inputFile.tellg();
    inputFile.seekg(0, std::ios::beg);

    // 确定文件中 double 的数量
    std::size_t numDoubles = fileSize / sizeof(double);

    // 创建一个与文件中 double 数量相同的向量
    std::vector<double> buffer(numDoubles);

    // 读取文件内容到缓冲区
    if (inputFile.read(reinterpret_cast<char*>(buffer.data()), fileSize)) {
        // std::cout << "文件读取成功" << std::endl;

        // 处理读取的数据（例如，输出前10个 double）
        // for (std::size_t i = 0; i < 10 && i < buffer.size(); ++i) {
        //     std::cout << buffer[i] << " ";
        // }
        std::cout << std::endl;
    } else {
        std::cerr << "读取文件时发生错误" << std::endl;
    }

    // 关闭文件
    inputFile.close();
    return buffer;
}

int compute_jacobi_set(
    const std::vector<double>& f, 
    const std::vector<double>& g, 
    const std::vector<std::pair<int, int> >& edges, 
    std::vector<std::vector<int> >& adjacency, 
    std::map<std::pair<int, int>, int>& is_critical, // 通过引用传递
    std::map<std::pair<int, int>, int>& decp_is_critical,
    double epsilon = 0,

    int type = 0) 
{

    std::map<std::pair<int, int>, int>* temp_is_critical;

    if(type==0){
        temp_is_critical = &is_critical;
    }
    else{
        temp_is_critical = &decp_is_critical;
    }
    for (int i =0;i<edges.size();i++) {
        std::pair<int, int> edge = edges[i];
        int a = edge.first;
        int b = edge.second;
        
        if (std::abs(g[b] - g[a]) == 0) {
            (*temp_is_critical)[edge] = 0;
            continue;
        }

        double h_lambda = (f[a] - f[b]) / (g[b] - g[a]);

        std::sort(adjacency[a].begin(), adjacency[a].end());
        std::sort(adjacency[b].begin(), adjacency[b].end());

        // 用于存储交集结果的 vector
        std::vector<int> v;

        // 使用 std::set_intersection 找到重复部分
        std::set_intersection(
            adjacency[a].begin(), adjacency[a].end(), 
            adjacency[b].begin(), adjacency[b].end(), 
            std::back_inserter(v)
        );
        // std::cout<<v.size()<<std::endl;
        if (v.size() != 2) continue;

        int v1 = v[0];
        int v2 = v[1];
        
        if ((f[v1] + h_lambda * g[v1] > f[a] + h_lambda * g[a] && f[v2] + h_lambda * g[v2] > f[a] + h_lambda * g[a]) || 
            (f[v1] + h_lambda * g[v1] < f[a] + h_lambda * g[a] && f[v2] + h_lambda * g[v2] < f[a] + h_lambda * g[a])) 
        {
            
            (*temp_is_critical)[edge] = 1;
            
            if(type==1 && is_critical[edge] == 0){
                false_edges.insert(edge);
            }
            
        } else {
            (*temp_is_critical)[edge] = 0;
            if(type==1 && is_critical[edge] == 1){
                false_edges.insert(edge);
            }
        }
    }

    return 0;
}

int compute_jacobi_set_local(std::pair<int, int> edge) {
    // 检查 (b, a) 是否存在于 edges 中
    int a = edge.first;
    int b = edge.second;

    // 检查 decp_g[b] 是否等于 decp_g[a]
    if (decp_g[a] == decp_g[b]) {
        return 0;
    } else {
        double h_lambda = (decp_f[a] - decp_f[b]) / (decp_g[b] - decp_g[a]);

        // 计算 a 和 b 的邻接点的交集
        std::vector<int> intersection;
        std::set_intersection(adjacency[a].begin(), adjacency[a].end(),
                              adjacency[b].begin(), adjacency[b].end(),
                              std::back_inserter(intersection));

        if (intersection.size() != 2) {
            return 0;
        }

        int v1 = intersection[0];
        int v2 = intersection[1];

        if ((decp_f[v1] + h_lambda * decp_g[v1] > decp_f[a] + h_lambda * decp_g[a] &&
             decp_f[v2] + h_lambda * decp_g[v2] > decp_f[a] + h_lambda * decp_g[a]) ||
            (decp_f[v1] + h_lambda * decp_g[v1] < decp_f[a] + h_lambda * decp_g[a] &&
             decp_f[v2] + h_lambda * decp_g[v2] < decp_f[a] + h_lambda * decp_g[a])) {
            return 1;
        } else {
            return 0;
        }
    }
}
int main(){
    input_f = getdata("f.bin");
    input_g = getdata("g.bin");
    decp_f = getdata("noisy_f.bin");
    decp_g = getdata("noisy_g.bin");
    double bound = 0.001;
    int width = 100;
    int height = 100;
    std::cout<<input_f.size()<<std::endl;


    std::vector<std::pair<int, int> > edges;
    

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            std::vector<int> temp;
            if (0 <= i - 1 && i-1< height && 0 <=(i - 1) * width + j && (i - 1) * width + j<width*height){
                temp.push_back((i - 1) * width + j);
                
                if ((std::find(edges.begin(), edges.end(), std::make_pair((i - 1) * width + j, i * width + j)) == edges.end()) &&
                    (std::find(edges.begin(), edges.end(), std::make_pair(i * width + j, (i - 1) * width + j)) ==edges.end()) ) {
                    edges.push_back(std::make_pair((i - 1) * width + j, i * width + j));
                }
            }

            if (0<= i - 1 && i-1 < height   && j + 1 < width && 0<=(i - 1) * width + j + 1 && (i - 1) * width + j + 1<width*height) {
                temp.push_back((i - 1) * width + j + 1);
                if((i - 1) * width + j + 1<0){
                    std::cout<<"chcuo2"<<std::endl;
                }
                if ((std::find(edges.begin(), edges.end(), std::make_pair((i - 1) * width + j + 1, i * width + j)) ==edges.end()) &&
                    (std::find(edges.begin(), edges.end(), std::make_pair(i * width + j, (i - 1) * width + j + 1)) ==edges.end())) {
                    edges.push_back(std::make_pair((i - 1) * width + j + 1, i * width + j));
                }
            }
            if (j - 1 >= 0  && 0<=i * width + j - 1 && i * width + j - 1<width*height) {
                temp.push_back(i * width + j - 1);
                
                if ((std::find(edges.begin(), edges.end(), std::make_pair(i * width + j - 1, i * width + j)) ==edges.end()) &&
                    (std::find(edges.begin(), edges.end(), std::make_pair(i * width + j, i * width + j - 1)) ==edges.end())) {
                    edges.push_back(std::make_pair(i * width + j - 1, i * width + j));
                }
            }

            if (j + 1 < width && 0<=i * width + j + 1 && i * width + j + 1<width*height) {
                
                temp.push_back(i * width + j + 1);
               
                if ((std::find(edges.begin(), edges.end(), std::make_pair(i * width + j + 1, i * width + j)) ==edges.end())&&
                     (std::find(edges.begin(), edges.end(), std::make_pair(i * width + j, i * width + j + 1)) ==edges.end())) {
                    edges.push_back(std::make_pair(i * width + j + 1, i * width + j));
                }
            }

            if (i + 1 < height && j - 1 >= 0 && 0<=(i + 1) * width + j - 1 && (i + 1) * width + j -1 < width*height) {
                temp.push_back((i + 1) * width + j - 1);
                 
                if ((std::find(edges.begin(), edges.end(), std::make_pair((i + 1) * width + j - 1, i * width + j)) ==edges.end()) &&
                    (std::find(edges.begin(), edges.end(), std::make_pair(i * width + j, (i + 1) * width + j - 1)) ==edges.end())) {
                    edges.push_back(std::make_pair((i + 1) * width + j - 1, i * width + j));
                }
            }

            if (i + 1 < height && 0<=(i + 1) * width + j && (i + 1) * width + j<width*height) {
                temp.push_back((i + 1) * width + j);
               
                if ((std::find(edges.begin(), edges.end(), std::make_pair((i + 1) * width + j, i * width + j)) ==edges.end()) &&
                    (std::find(edges.begin(), edges.end(), std::make_pair(i * width + j, (i + 1) * width + j)) ==edges.end())) {
                    edges.push_back(std::make_pair((i + 1) * width + j, i * width + j));
                }
            }

            adjacency.push_back(temp);


    }}
    
    
    compute_jacobi_set(input_f, input_g, edges, adjacency, is_critical, decp_is_critical, 0, 0);
    compute_jacobi_set(decp_f, decp_g, edges, adjacency, is_critical, decp_is_critical,0, 1);
    // for(int i=0;i<edges.size();i++){
    //     std::cout<<is_critical[i]<<", "<<decp_is_critical[i]<<std::endl;
    // }
    // std::cout<<false_edges.size()<<std::endl;
    
    // 统计每个元素出现的次数
    for (const auto& edge : false_edges) {
        frequency_map[edge.first]++;
        frequency_map[edge.second]++;
    }

    // 将元素及其出现次数存储到 std::vector<std::pair<int, int>> 中
    std::vector<std::pair<int, int> > frequency_vector(frequency_map.begin(), frequency_map.end());

    // 按照出现次数对 std::vector 进行排序（降序）
    std::sort(frequency_vector.begin(), frequency_vector.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return a.second > b.second;
    });
    // std::cout << "元素出现次数排序结果: " << std::endl;
    // for (const auto& elem : frequency_vector) {
    //     std::cout << "元素: " << elem.first << " 出现次数: " << elem.second << std::endl;
    // }
    while(false_edges.size()>0){
        frequency_map.clear();
        
        for (const auto& edge : false_edges) {
            // std::cout<<edge.first<<std::endl;
            frequency_map[edge.first]++;
            frequency_map[edge.second]++;
        }
        
        // 将元素及其出现次数存储到 std::vector<std::pair<int, int>> 中
        std::vector<std::pair<int, int> > frequency_vector(frequency_map.begin(), frequency_map.end());

        // 按照出现次数对 std::vector 进行排序（降序）
        std::sort(frequency_vector.begin(), frequency_vector.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            return a.second > b.second;
        });
        
        std::cout<<frequency_vector.size()<<std::endl;
        for(const auto& elem : frequency_vector){
            int v0 = elem.first;
            // std::cout<<v0<<std::endl;
            for(const int& i:adjacency[v0]){
                std::pair<int, int> edge;
                if(std::find(edges.begin(), edges.end(), std::make_pair(i, v0)) != edges.end()){
                    edge = std::make_pair(i, v0);
                }
                else{
                    edge = std::make_pair(v0, i);
                }
                int true_type = is_critical[edge];
                int false_type = compute_jacobi_set_local(edge);
                // std::cout<<true_type <<", "<<false_type<<std::endl;
                if(false_type!=true_type){
                    // std::cout<<"wrong"<<std::endl;
                    int a = edge.first;
                    int b = edge.second;
                    double noisy_lambda = (decp_f[a] - decp_f[b]) / (decp_g[b] -decp_g[a]);
                    double h_lambda = (input_f[a] - input_f[b]) / (input_g[b] -input_g[a]);
                    std::vector<int> v;
                    std::set_intersection(
                        adjacency[a].begin(), adjacency[a].end(), 
                        adjacency[b].begin(), adjacency[b].end(), 
                        std::back_inserter(v)
                    );
                    
                    if (v.size() != 2) continue;

                    int v1 = v[0];
                    int v2 = v[1];
                    int sign_v1, sign_v2, decp_sign_v1, decp_sign_v2;
                    sign_v1 = (input_f[v1] + h_lambda * input_g[v1] - (input_f[a] + h_lambda * input_g[a]))<0?1:0;
                    sign_v2 = (input_f[v2] + h_lambda * input_g[v2] - (input_f[a] + h_lambda * input_g[a]))<0?1:0;
                    decp_sign_v1 = (decp_f[v1] + noisy_lambda * decp_g[v1] - (decp_f[a] + noisy_lambda * decp_g[a]))<0?1:0;
                    decp_sign_v2 = (decp_f[v2] + noisy_lambda * decp_g[v2] - (decp_f[a] + noisy_lambda * decp_g[a]))<0?1:0;
                    // check which vertex is wrong.
                    // if v1 is wrong: 
                    
                    if(sign_v1!=decp_sign_v1){
                        // std::cout<<"v1 is wrong"<<std::endl;
                        double threshold_f = decp_f[a] + noisy_lambda * decp_g[a] - noisy_lambda * decp_g[v1];
                        double threshold_g = (decp_f[a] + noisy_lambda * decp_g[a] - decp_f[v1]) / (noisy_lambda);
                        
                        if(input_f[v1] - bound<=threshold_f && threshold_f<=input_f[v1] + bound){
                        
                            if(threshold_f<decp_f[v1]){ 
                                
                                decp_f[v1] = (threshold_f + input_f[v1] -bound)/2;
                            }
                            else{ 
                                decp_f[v1] = (threshold_f + input_f[v1] + bound)/2;
                                }
                            }
                        else if(input_g[v1] - bound<=threshold_g<=input_g[v1] + bound){
                        // # noisy_is_critical[edge] = is_critical[edge]
                        // std::cout<<"这里改1"<<std::endl;
                        if(threshold_g<decp_g[v1]){
                            
                            decp_g[v1] = (threshold_g + input_g[v1] -bound)/2;}
                        else{
                            decp_g[v1] = (threshold_g + input_g[v1] + bound)/2;
                            }
                        }
                    }

                    else if(sign_v2!=decp_sign_v2){
                        // std::cout<<"v2 is wrong"<<std::endl;
                        double threshold_f = decp_f[a] + noisy_lambda * decp_g[a] - noisy_lambda * decp_g[v2];
                        double threshold_g = (decp_f[a] + noisy_lambda * decp_g[a] - decp_f[v2]) / (noisy_lambda);
                        if(input_f[v2] - bound<=threshold_f<=input_f[v2] + bound){
                        //  noisy_is_critical[edge] = is_critical[edge]
                        if(threshold_f<decp_f[v2]){ 
                            decp_f[v2] = (threshold_f + input_f[v2] -bound)/2;
                            }
                        else{ 
                            decp_f[v2] = (threshold_f + input_f[v2] + bound)/2;
                            }
                        }

                        else if(input_g[v2] - bound<=threshold_g<=input_g[v2] + bound){
                        // # noisy_is_critical[edge] = is_critical[edge]
                        if(threshold_g<decp_g[v2]){
                            decp_g[v2] = (threshold_g + input_g[v2] -bound)/2;}
                        else{
                            decp_g[v2] = (threshold_g + input_g[v2] + bound)/2;
                            }
                        }
                        
                    }
                         
            }
        }
        }
    false_edges.clear();
    compute_jacobi_set(decp_f, decp_g, edges, adjacency, is_critical, decp_is_critical,0, 1);
    
    
    }
    false_edges.clear();
    compute_jacobi_set(decp_f, decp_g, edges, adjacency, is_critical, decp_is_critical,0, 1);
    
    return 0;
}