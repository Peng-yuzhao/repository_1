#include <stdio.h>
#include <stdlib.h>

/* run this program using the console pauser or add your own getch, system("pause") or input loop */





void input_linear(int line,int row,float* dicision,float* constant,float** coe){
	int i,j,k,m,capacity;
	char *c,*t;
	t=(char*)malloc(sizeof(char)*33);
	printf("������������߱���ǰ��ϵ��(��%d��)\n",row);
	for(i=0;i<row;++i){
		scanf("%f",&dicision[i]); 
	}
	getchar();
	printf("������Լ��������%d*%dϵ������\n",line,row);
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
	printf("����������Լ�������ĳ�����(��%d��)\n",line);
	for(i=0;i<line;++i){
		scanf("%f",&constant[i]);
	}
	printf("������ɣ�\n");
}




void simplex_start(float* dicision,float* constant,float** coe,int* base_variable,int* else_variable,int* m,int length_dic,int length_con){
	            //base_variable��else_variable��������ļ�������0��ʼ��ʾ 
	int i,j,l;
	float k;
    float max=-1;
    for(i=0;i<length_dic-length_con;++i)
    	if(dicision[i]>max){
    	max=dicision[i];
		j=i;
		}
    if(max<0)   m[0]=-1;
    else{
    if(max==0)   m[0]=-2,m[1]=j;
	else{
	for(i=0;i<length_con;++i)
	   if(coe[i][j]!=0){
	   k=constant[i]/coe[i][j];
	   break;
	}
	if(i==length_con)   m[0]=-3;
	else{
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
}
}
}




void simplex_transform(float* dicision,float* constant,float** coe,int length_dic,int length_con){
	int i;
	for(i=0;i<length_con;++i)
	coe[i][length_dic]=constant[i];
	for(i=0;i<length_dic;++i)
	coe[length_con][i]=dicision[i];
}



void simplex_calculation(float** form,int* base_variable,int* else_variable,int m[2],int length_form,int length_con){
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
             if(max<-0.00001)   m[0]=-1;
             else{
			 if(max>=-0.00001&&max<=0.00001)    m[0]=-2,m[1]=j;
			 else{
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
	     i=base_variable[l];
     	 base_variable[l]=else_variable[j];
	     else_variable[j]=i;
	     m[0]=l,m[1]=j;      
}
}
}


 
 
int main(int argc, char *argv[]) {
	int a,b,m[2],i,j;
	printf("�������������Թ滮��׼�ͷ���������δ֪����������������Ϊmax��\n");
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
	if(m[0]==-1){
		printf("ǰ%d�������������зǻ�����ȡֵΪ0ʱ���õ����Ž�0\n",b-a);
		return 0;
	}
		if(m[0]==-2){
		i=m[1]+1;
		printf("����%d�������⣬ǰ%d��������������%d���ǻ������⣬���зǻ�����ȡֵΪ0ʱ���õ����Ž�0.��%d����������%d���ǻ�������ȡ����Լ������������ֵ\n",i,b-a,i,i,i);
		return 0;
	}
	    if(m[0]==-3){
	    i=m[1]+1;
	    printf("��%d������������%d���ǻ�������ȡ���ڵ���0������ֵ������ȡֵ�����������ʱ���õ����Ž������������\n",i,i);	
	    return 0;
		}
    simplex_transform(dicision,constant,coe,b,a);
    do   simplex_calculation(coe,base_variable,else_variable,m,b+1,a);
    while(m[0]>=0);
    if(m[0]==-1){
    	for(i=0;i<a;++i){
    	variable[base_variable[i]][0]=0;  
    	variable[base_variable[i]][1]=i;
		}                                     //0�����������1����ǻ����� 
    	for(i=0;i<b-a;++i)
    	variable[else_variable[i]][0]=1;
	printf("��������ȡֵΪ(");
		float base_number[a];
		j=0;
		if(variable[0][0]==0){
		base_number[j++]=coe[variable[0][1]][b];
		printf("%f",base_number[0]);
	}
	else   printf("0");
	for(i=1;i<b;++i){
		if(variable[i][0]==0){
			base_number[j++]=coe[variable[i][1]][b];
			printf(",%f",base_number[j-1]);
		}
		else   printf(",0");
	}
    for(i=0,j=0,x=0;i<b;++i)
    	if(variable[i][0]==0)
    		x=x+base_number[j++]*dicision[base_variable[variable[i][1]]];
	printf(")ʱ��ȡ�����Ž�%f",x);
	return 0;
    }
    if(m[0]==-2){
    	for(i=0;i<a;++i){
    	variable[base_variable[i]][0]=0;  
    	variable[base_variable[i]][1]=i;
		}                                                //0�����������1����ǻ����� 
    	for(i=0;i<b-a;++i)
    	variable[else_variable[i]][0]=1;
	    printf("��%d������Ҳ������Ϊ���������ʸ����������������Ž⣬�ҵ���������ȡֵΪ(",else_variable[m[1]]+1);
	    j=0;
	    float base_number[a];
	if(variable[0][0]==0){
		base_number[j++]=coe[variable[0][1]][b];
		printf("%f",base_number[0]);
	}
	else   printf("0");
	for(i=1;i<b;++i){
		if(variable[i][0]==0){
			base_number[j++]=coe[variable[i][1]][b];
			printf(",%f",base_number[j-1]);
		}
		else   printf(",0");
	}
    for(i=0,j=0,x=0;i<b;++i)
    	if(variable[i][0]==0)
    		x=x+base_number[j++]*dicision[base_variable[variable[i][1]]];
	printf(")ʱ��ȡ��һ�����Ž�%f",x);
	return 0;
}
	return 1;
}
