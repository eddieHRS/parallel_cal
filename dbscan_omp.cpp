#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>

#include <algorithm>
#include <string>
#include <time.h>
#include <stdlib.h>  
#include <math.h>
using namespace std;
#include <string.h>  
#include <string>
#include <omp.h>  
using std::string; 


class point{  
public:  
    int dim;
	vector<double> value;
	int cluster;  
    int pointType;//1 noise 2 border 3 core  
    int visited;  
	point (){
		value.clear();
		cluster=-1;
		pointType=-1;
		visited=0;
	}
	double distance(point *p){
		double dist=0;
		for(int i=0;i<dim;i++)
			dist=dist+(value[i]-p->value[i])*(value[i]-p->value[i]);
		dist=sqrt(dist);
		return dist;
	}
	double caldis(point* p, int left,int right){
		double dist = 0;
		for(int i=left;i<=right;i++)
			dist=dist+(value[i]-p->value[i])*(value[i]-p->value[i]);
		return dist;
	}
};  

int thread_count = 16;
vector<point *> points; 
//vector<point *> &points = init_points;
int global_cluster=0;
const int Eps=10;
const int Minpts=200;

void SplitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
{
  std::string::size_type pos1, pos2;
  pos2 = s.find(c);
  pos1 = 0;
  while(std::string::npos != pos2)
  {
    v.push_back(s.substr(pos1, pos2-pos1));
 
    pos1 = pos2 + c.size();
    pos2 = s.find(c, pos1);
  }
  if(pos1 != s.length())
    v.push_back(s.substr(pos1));
}

void readFile(string filename){
	points.clear();
	ifstream inf;
	inf.open(filename.c_str()); 
	string s;
	while (getline(inf, s))
    {
		vector<string> a;
		SplitString(s,a,",");
		point *p=new point();
		p->dim=a.size()-1;
		for(int i=1;i<a.size();i++){
			double tmp;
			stringstream stream(a[i]);  
			stream>>tmp;  
			p->value.push_back(tmp);
		}
		points.push_back(p);
	}      
}

vector<point *> computeEps(point *p){
	vector<point *> Neighbor;
	Neighbor.clear();
#pragma omp parallel for shared(p,points)
	for(int i=0;i<points.size();i++){
		//对于那些已经确定集群的点不进行计算 global_cluster
		if(points[i] -> cluster == -1){
			double d1,d2,d3,d4,dis_all;
			#pragma private(d1,d2,d3,d4,dis_all)
			d1 = p->caldis(points[i],0,7);
			if(sqrt(d1) <= Eps){
				d2 = p->caldis(points[i],8,15);
				if(sqrt(d1 + d2) <= Eps){
					d3 = p->caldis(points[i],16,23);
					if(sqrt(d1+d2+d3) <= Eps){
						d4 = p->caldis(points[i],24,31);
						dis_all = sqrt(d1+d2+d3+d4);
						if(dis_all!= 0 && dis_all<=Eps){
							#pragma omp critical
							Neighbor.push_back(points[i]);
						}
					}
				}
			}
		}		
	}
	return Neighbor;
}

void expendSeed(vector<point*>& Seed){
//#pragma omp parallel for num_threads(8)
	for(int index = 0; index < Seed.size();index++){
		point *q=Seed[index];
		if(q->cluster == -1){
			q->cluster=global_cluster;
			if(q->pointType == 1)
				continue;
			vector<point *> Neighbor=computeEps(q);
			if(Neighbor.size()>=Minpts){
				q->visited=1;
				q->pointType=3;
				expendSeed(Neighbor);
			}
			else{
			 	q->visited=1;
			 	q->pointType=2;
			}
		}
	}
}

void dbscan(){
//#pragma omp parallel for num_threads(thread_count)
	for(int i = 0; i < points.size();i++){
		point* p = points[i];
		if(p -> cluster == -1){
			vector<point *> Neighbor=computeEps(p);
			if(Neighbor.size() >= Minpts){
				p->visited=1;
				p->pointType=3;
				p->cluster=++global_cluster;
				expendSeed(Neighbor);
			}
			//noise piont
			else{
			 	p->visited=1;
				p->pointType=1;
			}
		}
	}
}

int main(){
	//double t1 = omp_get_wtime();
	readFile("input.txt");
	double start = omp_get_wtime();
	dbscan();
	double end = omp_get_wtime();
	printf("%lf\n",end - start);
	ofstream fout("output.txt");     
	for(int i=0;i<points.size();i++)
		fout <<points[i]->cluster<< endl;
	fout.close(); 
	//double t2 = omp_get_wtime();
	//printf("all time:%lf\n",t2-t1);
	// double s = omp_get_wtime();
	
	// for(int i=0;i<points.size();i++)
	// 	delete points[i];
	// double e = omp_get_wtime();
	// printf("e - s%lf\n",e - s);

	return 0;
}
