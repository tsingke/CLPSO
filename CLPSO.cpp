//CLPSO 10.09.28
#include<iostream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
using namespace std;
const double PI = acos(-1.0);
const double MAXF = pow(2.0,1023);//maxdouble
const int DMS = 30; // dimension
const int ps = 40; // population size
const int mg =6000; // max generation
const double Vmax = 100;
const double Vmin = -100;
const double MINF = -1000000000;
const int m =7;  //abcd
double w0=0.9; 
double w1=0.4;
double c1=1.49445,c2=1.49445;
double w;   // inertia weight
double gbest[DMS];
double bestfitness;
double Xup[DMS];
double Xdown[DMS]; // x的范围;
double Pc[ps];
int flag[ps];
int fid[DMS];
int f1,f2;
int k,i,d;// k generation counter ; i partical id; d dimension id;
ofstream out;
char name[25];
char outname[25];
char txt[]=".txt";
char txtt[]="\\";
int indexname;
int gindex;

struct partical{
	double x[DMS];
	double v[DMS];
	double pbest[DMS];
	double fitness;
	double pfitness;
	bool feasible; //是否越界
};
partical P[ps];
void setdata(){
	FILE * IN;
	IN=fopen("D://range//in1.txt","r");        //abcd
	for(int t=0;t<DMS;t++){
		fscanf(IN,"%lf %lf",&Xdown[t],&Xup[t]);
	}	
}

void fitness1(int id){//Sphere unimodal
	P[id].fitness = 0; // 根据测试函数修改；
	for(int t=0;t<DMS;t++){
		P[id].fitness+=P[id].x[t]*P[id].x[t];
	}
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =1/P[id].fitness ;
}
void fitness2(int id){     //f2          
	P[id].fitness = 0; // 根据测试函数修改；
	double temp=1;
	for(int t=0;t<DMS;t++){
		P[id].fitness+=P[id].x[t];
		temp*=abs(P[id].x[t]);
	}
	P[id].fitness +=temp;
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =1/P[id].fitness ;
}
void fitness3(int id){
	P[id].fitness = 0; // 根据测试函数修改；
	double temp;
	for(int t=0;t<DMS;t++){
		temp=0;
		for(int d=0;d<=t;d++){
		temp+=P[id].x[d]*P[id].x[d];
		}
		P[id].fitness =temp*temp;
	}
	
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =1/P[id].fitness ;
}
void fitness4(int id){
	P[id].fitness =0;
	double temp=-1;
	double tempv;
	for(int t=0;t<DMS;t++){
		tempv=abs(P[id].x[t]);
		if(temp<tempv){
			temp=tempv;
		}
	}
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =temp;
}
void fitness5(int id){//Rosenbrock  simple multimodal
	P[id].fitness = 0; // 根据测试函数修改；
	for(int t=0;t<DMS-1;t++){
		P[id].fitness+=100*pow((P[id].x[t+1]-pow(P[id].x[t],2.0)),2.0)+pow((P[id].x[t]-1),2.0);
	}	
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =1/P[id].fitness ;
}
void fitness6(int id){
	P[id].fitness = 0; // 根据测试函数修改；
	for(int t=0;t<DMS;t++){
		P[id].fitness+=pow(floor(P[id].x[t]+0.5),2.0);
	}
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =1/P[id].fitness ;
}
void fitness7(int id){
	P[id].fitness = 0; // 根据测试函数修改；
	for(int t=0;t<DMS;t++){
		P[id].fitness+=t*pow(P[id].x[t],4.0);
	}
	P[id].fitness+=(rand()%1000/1000.0);
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =1/P[id].fitness ;
}
void fitness8(int id){
	P[id].fitness = 0; // 根据测试函数修改；
	for(int t=0;t<DMS;t++){
		P[id].fitness+=(-P[id].x[t]*sin(sqrt(abs(P[id].x[t]))));
	}
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =1/P[id].fitness ;
}
void fitness9(int id){//Rastrigin unrotated multimodal
	P[id].fitness = 0; // 根据测试函数修改；
	for(int t=0;t<DMS;t++){
		P[id].fitness+=(pow(P[id].x[t],2.0)-10*cos(2*PI*P[id].x[t])+10);
	}
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =1/P[id].fitness ;
}
void fitness10(int id){//Ackley unrotated multimodal
	P[id].fitness = 0; // 根据测试函数修改；
	double temp0=0;
	double temp1=0;
	for(int t=0;t<DMS;t++){
		temp0+=pow(P[id].x[t],2.0);
		temp1+=cos(2*PI*P[id].x[t]);
	}
	P[id].fitness = (-20*(exp(-0.2*sqrt(temp0/DMS)))-exp(temp1/(double)DMS)+20+exp(1.0));
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =1/P[id].fitness ;
}
void fitness11(int id){//Griewanks unrotated multimodal
	P[id].fitness = 0; // 根据测试函数修改；
	double temp0=0;
	double temp1=1;
	for(int t=0;t<DMS;t++){
		temp0+=pow(P[id].x[t],2.0);
		temp1*=cos(P[id].x[t]/sqrt((double)t));
	}
	P[id].fitness =temp0/4000.0-temp1+1;
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =1/P[id].fitness ;
}
void fitness13(int id){
	P[id].fitness = 0; // 根据测试函数修改；
	double temp=0;
	double temp0=0;
	for(int t=0;t<DMS-1;t++){
		temp+=pow((P[id].x[t]-1),2.0)*(1+pow(sin(3*PI*P[id].x[t+1]),2.0));
	}	
	for(int t=0;t<DMS;t++){
		if(P[id].x[t]>5){
			temp0+=100*pow((P[id].x[t]-5),4.0);
		}else if(P[id].x[t]<-5){
			temp0+=100*pow((-P[id].x[t]-5),4.0);
		}else {
			temp0+=0;
		}
	}
	P[id].fitness = 0.1*(pow(sin(3*PI*P[id].x[0]),2.0)+temp
		+(P[id].x[DMS-1]-1)*(1+pow(sin(2*PI*P[id].x[DMS-1]),2.0)))+temp0;
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness =1/P[id].fitness ;
}
void Fitness5(int id){
	double a,b;
	a=0.5;
	b=3;
	int kmax=20;
	P[id].fitness =0;
	double temp=0;
	for(int t=0;t<DMS;t++){
		for(int r=0;r<=kmax;r++){
			P[id].fitness +=pow(a,k)*cos(2*PI*pow(b,k)*(P[id].x[t]+0.5));
		}
	}
	for(int r=0;r<=kmax;r++){
		temp+=pow(a,k)*cos(2*PI*pow(b,k)*0.5);
	}
	P[id].fitness =P[id].fitness -DMS*temp;

	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness = 1/P[id].fitness ;

}
void Fitness7(int id){ // noncontinuous Rastrigin's function;
	double y[DMS];
	P[id].fitness =0;
	for(int t=0;t<DMS;t++){
		if(abs(P[id].x[t])<0.5)y[t]=P[id].x[t];
		else y[t]=int(2*P[id].x[t]+0.5)/2.0;
		P[id].fitness +=(y[t]*y[t]-10*cos(2*PI*y[t])+10);
	}
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness = 1/P[id].fitness ;

}
void Fitness8(int id){
//Schwefel's function;   
	P[id].fitness =0;
	for(int t=0;t<DMS;t++){
		P[id].fitness +=P[id].x[t]*sin(sqrt(abs(P[id].x[t])));
	}
	P[id].fitness=418.9829*DMS-P[id].fitness ;
	if(P[id].fitness ==0)P[id].fitness =MAXF;
	else 
	P[id].fitness=1/P[id].fitness;
}
void initial(){
	int temp;
	double tempv=MINF;
	
	for(int i=0;i<ps;i++){
		P[i].feasible =true;
		for(int j=0;j<DMS;j++){
			P[i].v[j]=((rand()%1000)/1000.0 )*(Vmax-Vmin)+Vmin;
			P[i].x[j]=P[i].pbest [j]=((rand()%1000)/1000.0 )*(Xup[j]-Xdown[j])+Xdown[j];
		}
		fitness1(i);      //abcd
		P[i].pfitness =P[i].fitness ;
		if(tempv<P[i].pfitness){
			tempv=P[i].pfitness;
			temp=i;
		}
		Pc[i]=0.05+0.45*(exp((10.0*i)/(ps-1))-1)/(exp(10.0)-1);
		flag[i]=0;
	}
	bestfitness=tempv;
	gindex=temp;
	for(int j=0;j<DMS;j++){
		gbest[j]=P[temp].x[j] ;
	}

}
void Learning(int id,int f,int did){
	P[id].v [did]=w*P[id].v [did]+c1*((rand()%1000)/1000.0)*(P[f].pbest[did]-P[id].x[did]);
	if(P[id].v [did]<Vmin)P[id].v [did]=Vmin;
	if(P[id].v [did]>Vmax)P[id].v [did]=Vmax;
	P[id].x[did]=P[id].x[did]+P[id].v [did];
	if(P[id].x[did]<Xdown[did]){
		P[id].x[did]=Xdown[did];
		P[id].feasible =false;
	}
	if(P[id].x[did]>Xup[did]){
		P[id].x[did]=Xup[did];
		P[id].feasible =false;
	}
}
void learning(int id,int did){
	P[id].v [did]=w*P[id].v [did]+c2*((rand()%1000)/1000.0)*(P[gindex].pbest [did]-P[id].x[did])
		+c1*((rand()%1000)/1000.0)*(P[id].pbest[did]-P[id].x[did]);
		
	if(P[id].v [did]<Vmin)P[id].v [did]=Vmin;
	if(P[id].v [did]>Vmax)P[id].v [did]=Vmax;
	P[id].x[did]=P[id].x[did]+P[id].v [did];
	if(P[id].x[did]<Xdown[did]){
		P[id].x[did]=Xdown[did];
		P[id].feasible =false;
	}
	if(P[id].x[did]>Xup[did]){
		P[id].x[did]=Xup[did];
		P[id].feasible =false;
	}
}
void print(){
	printf("%.50lf \n",1/bestfitness);
	out<<1/bestfitness<<endl;;
}
void updatefitness(int id){
	if(P[id].fitness >P[id].pfitness){
		for(int t=0;t<DMS;t++){
			P[id].pbest [t]=P[id].x[t];
		}
		flag[id]=0;
		P[id].pfitness =P[id].fitness ;
		if(P[id].fitness >bestfitness){
			for(int t=0;t<DMS;t++){
				gbest[t]=P[id].x[t];
			}
			bestfitness=P[id].fitness ;
		}
	}else{
		flag[id]++;
	}
}
void Update(int id){
	if(P[id].fitness >P[id].pfitness){
		for(int t=0;t<DMS;t++){
			P[id].pbest [t]=P[id].x[t];
		}
		P[id].pfitness =P[id].fitness ;
		if(P[id].fitness >bestfitness){
			gindex=id;
			for(int t=0;t<DMS;t++){
				gbest[t]=P[id].x[t];
			}
			bestfitness=P[id].fitness ;
		}
	}
}
int main(){
	//clock_t start,end;
	srand(time(NULL));
	setdata();
	int test[3];
	
	for(int ddd=1;ddd<=30;ddd++){ // 运行多少次times
//	start=clock();
	//PSO
		test[0]=test[1]=test[2]=0;
		out.clear ();
		out.close ();
		strcpy(outname,"D://f1//");//abcd 关键位置。。。
		itoa(ddd,name,10);
		strcat(outname,name);
		strcat(outname,txt);
		out.open (outname);
		initial();
	
/*	for(k=0;k<mg;k++){
	//	w=(w0-w1)*k/mg+w1;
		w=(w0-w1)*pow(((double)k/mg-1),2.0)+w1;
	//	w=w1;
		
		for(i=0;i<ps;i++){
			for(d=0;d<DMS;d++){
				learning(i,d);			
			}
		//	if(P[i].feasible ){
			fitness1(i);//abcd
			Update(i);
			test[2]++;
		//	}
		}
		print();
	}*/
	//CLPSO
	///*
	for(k=0;k<mg;k++){
		if(test[2]>200000)break;   //abcd
	//	w=(w0-w1)*pow(((double)k/mg-1),2.0)+w1;
		w=w0*(w0-w1)*k/mg;
	//	w=0;
		for(i=0;i<ps;i++){
			if(flag[i]>=m){//一直没有更新到pbest的话，则选择另外的学习方式；
				for(d=0;d<DMS;d++){
					if((rand()%1000/1000.0)<Pc[i]){
						f1=floor((rand()%1000/1000.0)*ps);
						f2=floor((rand()%1000/1000.0)*ps);
						if(P[f1].pfitness>P[f2].pfitness){
							fid[d]=f1;
						}else{
							fid[d]=f2;
						}
					}else{
						fid[d]=i;
					}				
				}
				flag[i]=0;
				for(d=0;d<DMS;d++){
					Learning(i,fid[d],d);				
				}
				test[0]++;
			}else{//优先学习自己的最优解
				int temp=rand()%DMS;
				for(d=0;d<DMS;d++){
					if(d==temp)fid[d]=rand()%ps;
					else fid[d]=i;
				}
				for(d=0;d<DMS;d++){
					Learning(i,fid[d],d);			
				}
				test[1]++;
			}

		//    if(P[i].feasible ){
				fitness1(i); //abcd;
				updatefitness(i);
				test[2]++;
		//	}else P[i].feasible=true; 
		}
		print();
		
	}//*/
//	end=clock();
//	cout<<end-start<<endl;	
	
	cout<<test[2]<<endl;
	}
	cout<<(double)test[0]/(double)test[1]<<endl;


	cout<<endl;

	return 0;
}