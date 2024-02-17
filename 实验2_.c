#include <stdio.h>
#include <stdlib.h>

/* run this program using the console pauser or add your own getch, system("pause") or input loop */





void input_linear(int line,int row,float* dicision,float* constant,float** coe){
	int i,j,k,m,capacity;
	char *c,*t;
	t=(char*)malloc(sizeof(char)*33);
	printf("请依次输入决策变量前的系数(共%d个)\n",row);
	for(i=0;i<row;++i){
		scanf("%f",&dicision[i]); 
	}
	getchar();
	printf("请输入约束条件的%d*%d系数矩阵\n",line,row);
	for(i=0;i<line;++i){
		capacity=10*row;
		c=(char*)malloc(sizeof(char)*capacity);
		for(m=0;;++m){
			if(m==capacity-1){
			capacity=capacity+5*row;
			c=(char*)realloc(c,sizeof(char)*capacity);
		}
				if((c[m]=getchar())=='\n')   break;
		}
		for(j=0,m=0;j<row;++j){
		    for(k=0;c[m]!=' ';++k,++m){
		    if(c[m]=='\n')   break;
		    t[k]=c[m];
		    }
		    t=(char*)realloc(t,sizeof(char)*(k+1));
		    coe[i][j]=atof(t);
		    if(c[m]=='\n')   break;
		    do   ++m;
		    while(c[m]==' ');
		    if(c[m]=='\n')   break;
		}
	}
	printf("请依次输入约束条件的常数项(共%d个)\n",line);
	for(i=0;i<line;++i){
		scanf("%f",&constant[i]);
	}
	printf("输入完成！\n");
}




void simplex_start(float* dicision,float* constant,float** coe,int* base_variable,int* else_variable,int* m,int length_dic,int length_con){
	            //base_variable、else_variable按照数组的计数，从0开始表示 
	int i,j,l;
	float k;
    float max=-1;
    for(i=0;i<length_dic-length_con;++i)
    	if(dicision[i]>max){
    	max=dicision[i];
		j=i;
		}
    if(max<0){
    	printf("前%d个变量，即所有非基变量取值为0时，取得最优解0\n",length_dic-length_con);
		exit(0);
	}
    if(max==0){
    	printf("除第%d个变量外，前%d个变量，即除第%d个非基变量外，所有非基变量取值为0时，取得最优解0.第%d个变量，即%d个非基变量可取满足约束条件的任意值\n",j+1,length_dic-length_con,j+1,j+1,j+1);
		exit(0);
	}

	for(i=0;i<length_con;++i)
	   if(coe[i][j]!=0){
	   k=constant[i]/coe[i][j];
	   break;
	}
	if(i==length_con){
		printf("第%d个变量，即第%d个非基变量可取大于等于0的任意值。当其取值趋近于无穷大时，取得最优解趋近于无穷大\n",j+1,j+1);	
	    exit(0);
	}
	for(l=i++;i<length_con;++i)
		if(coe[i][j]!=0&&constant[i]/coe[i][j]<k){
			k=constant[i]/coe[i][j];
			l=i;
		}
	i=base_variable[l];
	base_variable[l]=else_variable[j];
	else_variable[j]=i;
	m[0]=l,m[1]=j;
}




void simplex_transform(float* dicision,float* constant,float** coe,int length_dic,int length_con){
	int i;
	for(i=0;i<length_con;++i)
	coe[i][length_dic]=constant[i];
	for(i=0;i<length_dic;++i)
	coe[length_con][i]=dicision[i];
}



int simplex_calculation(float** form,int* base_variable,int* else_variable,int m[2],int length_form,int length_con){
	 float max=-1,standard,k;
	 int i,j,l;
	 standard=form[m[0]][base_variable[m[0]]];
	 for(j=0;j<length_form;++j)
	 	 form[m[0]][j]=form[m[0]][j]/standard;
	 standard=form[length_con][base_variable[m[0]]];
     for(j=0;j<length_form-1;++j)
     	 form[length_con][j]=form[length_con][j]-form[m[0]][j]*standard;
	 for(i=0;i<length_con;++i){
	 	 if(i==m[0]||form[i][base_variable[m[0]]]==0)   continue;
	 standard=form[i][base_variable[m[0]]];
	     for(j=0;j<length_form;++j)
	     	 form[i][j]=form[i][j]-form[m[0]][j]*standard;
	     }
         for(i=0;i<length_form-length_con-1;++i)
             if(form[length_con][else_variable[i]]>max)
             j=i,max=form[length_con][else_variable[i]];
             if(max<-0.00001){
             	 m[0]=-1;
             	 return 0;
			 }
			 for(i=0;i<length_con;++i)
	         if(form[i][else_variable[j]]<-0.00001||form[i][else_variable[j]]>0.00001){
	         k=form[i][length_form-1]/form[i][else_variable[j]];
	         break;
             }
			 for(l=i++;i<length_con;++i)
	         if((form[i][else_variable[j]]<-0.00001||form[i][else_variable[j]]>0.00001)&&form[i][length_form-1]/form[i][else_variable[j]]<k){
			 k=form[i][length_form-1]/form[i][else_variable[j]];
			 l=i;
		}
			 if(max>=-0.00001&&max<=0.00001){
			 	  m[0]=-2,m[1]=j;
			 	  return l;
			 }
	     i=base_variable[l];
     	 base_variable[l]=else_variable[j];
	     else_variable[j]=i;
	     m[0]=l,m[1]=j;
		 return 0;      
}


void simplex_subsequence(float** form,float* dicision,int* base_variable,int* else_variable,int l,int m[2],int length_dic,int length_con){
	int i,j=0,variable[length_dic][2];
	float x,base_number[length_con];
    for(i=0;i<length_con;++i){
   	variable[base_variable[i]][0]=0;  
   	variable[base_variable[i]][1]=i;
	}                                              //0代表基变量，1代表非基变量 
    for(i=0;i<length_dic-length_con;++i)
   	variable[else_variable[i]][0]=1;
	if(m[0]==-1){
	printf("各个变量取值为(");
		if(variable[0][0]==0){
		base_number[j++]=form[variable[0][1]][length_dic];
		printf("%f",base_number[0]);
	}
	else   printf("0");
	for(i=1;i<length_dic;++i){
		if(variable[i][0]==0){
			base_number[j++]=form[variable[i][1]][length_dic];
			printf(",%f",base_number[j-1]);
		}
		else   printf(",0");
	}
	for(i=0,j=0,x=0;i<length_dic;++i)
    	if(variable[i][0]==0)
    		x=x+base_number[j++]*dicision[base_variable[variable[i][1]]];
	printf(")时，取得一个最优解%f",x);
	exit(0);
    }
    if(m[0]==-2){
	    printf("该问题有无穷多个最优解，且当各个变量取值为(");
	if(variable[0][0]==0){
		base_number[j++]=form[variable[0][1]][length_dic];
		printf("%f",base_number[0]);
	}
	else   printf("0");
	for(i=1;i<length_dic;++i){
		if(variable[i][0]==0){
			base_number[j++]=form[variable[i][1]][length_dic];
			printf(",%f",base_number[j-1]);
		}
		else   printf(",0");
	}
	for(i=0,j=0,x=0;i<length_dic;++i)
    	if(variable[i][0]==0)
    		x=x+base_number[j++]*dicision[base_variable[variable[i][1]]];
	printf(")时，取得一个最优解%f\n当第%d个变量取代第%d个变量作为基变量时，",x,else_variable[m[1]]+1,base_variable[l]+1);
	i=base_variable[l];
    base_variable[l]=else_variable[m[1]];
	else_variable[m[1]]=i;
	m[0]=l;
	simplex_calculation(form,base_variable,else_variable,m,length_dic+1,length_con);
	m[0]=-1;
	simplex_subsequence(form,dicision,base_variable,else_variable,0,m,length_dic,length_con);
	exit(0);
}
}
 
 
int main(int argc, char *argv[]) {
	int a,b,m[2],i;
	printf("请依次输入线性规划标准型方程数量和未知数数量（决策条件为max）\n");
	scanf("%d%d",&a,&b);
	float x,dicision[b],constant[a],**coe;
	coe=(float**)malloc(sizeof(float*)*(a+1));
	for(i=0;i<a;++i)
	coe[i]=(float*)malloc(sizeof(float)*(b+1));
	coe[a]=(float*)malloc(sizeof(float)*b);
	input_linear(a,b,dicision,constant,coe);
	int base_variable[a],else_variable[b-a],variable[b][2];
	for(i=0;i<b-a;++i)
	else_variable[i]=i;
	for(;i<b;++i)
	base_variable[i-b+a]=i;
	simplex_start(dicision,constant,coe,base_variable,else_variable,m,b,a);
    simplex_transform(dicision,constant,coe,b,a);
    do   
	i=simplex_calculation(coe,base_variable,else_variable,m,b+1,a);
    while(m[0]>=0);
    simplex_subsequence(coe,dicision,base_variable,else_variable,i,m,b,a);
}
