#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<omp.h>

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#define dynamic_count 64
int thread_count = 1;
double dist[1000010];
double rebuiltCoord[10010];
int kk[10];
//计算所有数据跟所有其他点之间的距离
//传入的数据包括 n=所有点的个数 coord=点的坐标
void Distance(const int n, const int dim, double* coord){
    //对所有的点遍历求距离
// #   pragma omp parallel for num_threads(thread_count) \
//         default(none) shared(n, dim, dist, coord) private(i, j, ki) schedule(static, 1)
#   pragma omp simd
    for(int i=0; i<n; i++)
        for(int j = i + 1; j < n; j ++){
            for(int ki = 0; ki < dim; ki ++) dist[i * n + j] += (coord[j*dim + ki] - coord[i*dim + ki]) * (coord[j*dim + ki] - coord[i*dim + ki]);
            dist[i * n + j] = sqrt(dist[i * n + j]);
        }
    
#   pragma omp simd
    for(int i = 0; i < n; i ++)
        for(int j = i - 1; j > -1; j --) dist[i * n + j] = dist[j * n + i];
    
}


double Chebyshev(const int n, const int k){
    double chebyshevSum = 0;

    int i, j, ki, ik, jk;
    double chebyshev;
#   pragma omp parallel for num_threads(thread_count) \
        default(none) shared(n, k, rebuiltCoord) private(i, j, ki, chebyshev, ik, jk) reduction(+:chebyshevSum) schedule(dynamic, dynamic_count)
        
    for(i=0; i<n; i++){
        ik = i * k;
        for(j= i + 1; j<n; j++){
            jk = j * k;
                // chebyshev = dis>chebyshev ? dis : chebyshev;
            chebyshev = fabs(rebuiltCoord[ik] - rebuiltCoord[jk]) + fabs(rebuiltCoord[ik+ 1] - rebuiltCoord[jk+ 1]);
            chebyshevSum = chebyshevSum + chebyshev;
        }  
    }
           
    

    
    return chebyshevSum;
}

// double fact(int n)
// {
// 	double factorial = 1;//局部变量声明
// 	for (; n >= 1; n--)
// 	{
// 		factorial *= n;
// 	}
// 	return factorial;//将计算结果返回
// }

double SumDistance(const int k, const int n, const int dim, double* coord, int* pivots){
    //n:点数 k:支撑点个数
    // double* rebuiltCoord = (double*)malloc(sizeof(double) * n * k);
    // int i;
    // for(i=0; i<n*k; i++){
    //     rebuiltCoord[i] = 0;
    // }


    // Rebuild coordinates. New coordinate of one point is its distance to each pivot.
    // 计算基于支撑点的坐标的位置
    int i, ki;
    //对所有点进行遍历
#   pragma omp parallel for num_threads(thread_count) \
        default(none) shared(n, k, rebuiltCoord, dist, pivots) private(i, ki) schedule(dynamic, dynamic_count)
    for(i=0; i<n; i++){
        //ki = 0, 1
        //新坐标的维数应该是k
        ki = 0;
        rebuiltCoord[i*k + ki] = dist[pivots[ki] * n + i] + dist[pivots[ki + 1] * n + i];
        ki = 1;
        rebuiltCoord[i*k + ki] = dist[pivots[ki - 1] * n + i] - dist[pivots[ki] * n + i];
    }
        
    

    // Calculate the sum of Chebyshev distance with rebuilt coordinates between every points
    // 计算重建坐标之后的点的切比雪夫距离和
    // double chebyshevSum = 0;
    // for(i=0; i<n; i++){
    //     int j;
    //     for(j=0; j<n; j++){
    //         double chebyshev = 0;
    //         int ki;
    //         for(ki=0; ki<k; ki++){
    //             double dis = fabs(rebuiltCoord[i*k + ki] - rebuiltCoord[j*k + ki]);
    //             chebyshev = dis>chebyshev ? dis : chebyshev;
    //         }
    //         chebyshevSum += chebyshev;
    //     }
    // }

    return Chebyshev(n, k);
}

void pivotSort(int k, int n, int *kk){
    int j, i, flag = 0;
    kk[k-1]++;
#   pragma omp simd
    for(j = k-1; j > 0; j --)
        if(unlikely(kk[j] > n - k + j )){
            kk[j - 1] = kk[j - 1] + 1;
            kk[j] = 0;
            flag = 1;
        } 
    
    if(unlikely(flag == 1))
        for(i = 1; i < k; i ++)
            if(kk[i] <= kk[i - 1]) kk[i] = kk[i - 1] + 1;
    
    
}

void SortSum(const int k, const int n, const int dim, double* coord, const int M, 
                 double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots){
    int i, j, kj, a, temp3, temp4, p = 0;
    double sum, temp1, temp2;
    // int* temp = (int*)malloc(sizeof(int) * k);
    Distance(n, dim, coord);
    // printf("原坐标距离矩阵：\n");
    // for(i = 0; i < n; i ++){
    //     for(int j = 0; j < n; j ++) printf("%lf ", dist[i * n + j]);
    //     printf("\n");
    // }

    //循环次数
    while(kk[0] < n - k + 1) {
        // printf("\n");
        // for(i = 0; i < k; i ++){
        //     printf("支撑点是%d ;\n", kk[i]);
        // }

        // for(i = 0; i < k; i ++) temp[i] = kk[i];
        sum = SumDistance(k, n, dim, coord, kk);

        if(likely(p > M || p == M)){
            // put data at the end of array
            maxDistanceSum[M] = sum;
            minDistanceSum[M] = sum;
        
            #   pragma omp simd
            for(kj=0; kj<k; kj++){
                maxDisSumPivots[M*k + kj] = kk[kj];
                minDisSumPivots[M*k + kj] = kk[kj];
            }
            

            // sort
        
            for(a=M; a>0; a--)
                if(maxDistanceSum[a] > maxDistanceSum[a-1]){
                    temp1 = maxDistanceSum[a];
                    temp2 = maxDistanceSum[a-1];
                    maxDistanceSum[a] = temp2;
                    maxDistanceSum[a-1] = temp1;
                    
                    for(kj=0; kj<k; kj++){
                        temp3 = maxDisSumPivots[a*k + kj];
                        temp4 = maxDisSumPivots[(a-1)*k + kj];
                        maxDisSumPivots[a*k + kj] = temp4;
                        maxDisSumPivots[(a-1)*k + kj] = temp3;
                    }
                }
                else break;
            
            for(a=M; a>0; a--)
                if(minDistanceSum[a] < minDistanceSum[a-1]){
                    temp1 = minDistanceSum[a];
                    temp2 = minDistanceSum[a-1];
                    minDistanceSum[a] = temp2;
                    minDistanceSum[a-1] = temp1;
                    
                    for(kj=0; kj<k; kj++){
                        temp3 = minDisSumPivots[a*k + kj];
                        temp4 = minDisSumPivots[(a-1)*k + kj];
                        minDisSumPivots[a*k + kj] = temp4;
                        minDisSumPivots[(a-1)*k + kj] = temp3;
                    }
                }
                else break;
            
        }
        else{
            maxDistanceSum[p] = sum;
            minDistanceSum[p] = sum;
            #   pragma omp simd
            for(kj=0; kj<k; kj++){
                maxDisSumPivots[p*k + kj] = kk[kj];
                minDisSumPivots[p*k + kj] = kk[kj];
            }
            
            
        }
        
        if(unlikely(kk[k - 1] == n - 1)){
            int kj;
            //支撑点
            for(kj=0; kj<k; kj++)
                printf("%d ", kk[kj]);
            
            putchar('\t');
            //最大距离支撑点坐标
            for(kj=0; kj<k; kj++)
                printf("%d ", maxDisSumPivots[kj]);
            
            //最大距离和
            printf("%lf\t", maxDistanceSum[0]);
        
            //最小距离和的支撑点
            for(kj=0; kj<k; kj++)
                printf("%d ", minDisSumPivots[kj]);
            
            //最小距离和
            printf("%lf\n", minDistanceSum[0]);
        }
        p ++;
        pivotSort(k, n, kk);
        
    } 
    
}


int main(int argc, char* argv[])
{
    // filename : input file namespace
    char* filename = (char*)"uniformvector-2dim-5h.txt";
    if( argc==2 ) {
        filename = argv[1];
    }  else if(argc != 1) {
        printf("Usage: ./pivot <filename>\n");
        return -1;
    }
    
    printf("请输入线程数：");
    scanf("%d", &thread_count);

    // M : number of combinations to store
    const int M = 1000;
    // dim : dimension of metric space 二维坐标(x, y)
    int dim;
    // n : number of points  点的个数
    int n;
    // k : number of pivots  支撑点个数
    int k;

    // Read parameter
    FILE* file = fopen(filename, "r");
    if( file == NULL ) {
        printf("%s file not found.\n", filename);
        return -1;
    }
    fscanf(file, "%d", &dim);
    fscanf(file, "%d", &n);
    fscanf(file, "%d", &k);
    printf("dim = %d, n = %d, k = %d\n", dim, n, k);

    // Start timing
    struct timeval start;

    // Read Data 读入点的坐标coord
    double* coord = (double*)malloc(sizeof(double) * dim * n);
    int i, ki, j;
#   pragma omp simd 
    for(i=0; i<n; i++)
         
        for(j=0; j<dim; j++)
            fscanf(file, "%lf", &coord[i*dim + j]);
        
    
    fclose(file);
    gettimeofday(&start, NULL);
    // maxDistanceSum : the largest M distance sum 
    double* maxDistanceSum = (double*)malloc(sizeof(double) * (M+1));
    // minDistanceSum : the smallest M distance sum
    double* minDistanceSum = (double*)malloc(sizeof(double) * (M+1));
    // minDisSumPivots : the bottom M pivots combinations
    int* minDisSumPivots = (int*)malloc(sizeof(int) * k * (M+1));
    // maxDisSumPivots : the top M pivots combinations
    int* maxDisSumPivots = (int*)malloc(sizeof(int) * k * (M+1));
#   pragma omp simd    
    for(i=0; i< M; i++){
        minDistanceSum[i] = __DBL_MAX__;
        maxDistanceSum[i] = 0;
         
        for(ki=0; ki<k; ki++){
            maxDisSumPivots[i*k + ki] = 0;
            minDisSumPivots[i*k + ki] = 0;
        }
    }
    
    // Distance(n, dim, coord);
    //初始化数组kk，使其从1开始递增
    for(i = k; i > 0; i --)
        for(j = k; j >= i; j --) kk[j] ++;
    

    SortSum(k, n, dim, coord, M, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);

    // End timing
    struct timeval end;
    gettimeofday (&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec-start.tv_sec)*1000.0+(end.tv_usec-start.tv_usec)/1000.0);

    // Store the result
    FILE* out = fopen("result.txt", "w");
   
    for(i=0; i<M; i++){
         
        for(ki=0; ki<k-1; ki++)
            fprintf(out, "%d ", maxDisSumPivots[i*k + ki]);
        
        fprintf(out, "%d\n", maxDisSumPivots[i*k + k-1]);
    }

    for(i=0; i<M; i++){
        
        for(ki=0; ki<k-1; ki++){
            fprintf(out, "%d ", minDisSumPivots[i*k + ki]);
        }
        fprintf(out, "%d\n", minDisSumPivots[i*k + k-1]);
    }
    fclose(out);

    // Log
    printf("max : ");
    for(ki=0; ki<k; ki++){
        printf("%d ", maxDisSumPivots[ki]);
    }
    printf("%lf\n", maxDistanceSum[0]);
    printf("min : ");
    for(ki=0; ki<k; ki++){
        printf("%d ", minDisSumPivots[ki]);
    }
    printf("%lf\n", minDistanceSum[0]);
    

    return 0;
    
    


}
