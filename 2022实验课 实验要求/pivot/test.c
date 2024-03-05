#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double dist[1000000];
double rebuiltCoord[10000];
//计算所有数据跟所有其他点之间的距离
//传入的数据包括 n=所有点的个数 coord=点的坐标
void Distance(const int n, const int dim, double* coord){
    int i, j, ki;
    //对所有的点遍历求距离
    for(i=0; i<n; i++){
        for(j = i + 1; j < n; j ++){
            for(ki = 0; ki < dim; ki ++) dist[i * n + j] += pow(coord[j*dim + ki] - coord[i*dim + ki] ,2);
            dist[i * n + j] = sqrt(dist[i * n + j]);
        }
    }
    for(i = 0; i < n; i ++){
        for(j = i - 1; j > -1; j --) dist[i * n + j] = dist[j * n + i];
    }
}

double Chebyshev(const int n, const int k){
    double chebyshevSum = 0;
    int i, j, ki;
    for(i=0; i<n; i++){
        for(j= i + 1; j<n; j++){
            double chebyshev = 0;
            for(ki=0; ki<k; ki++){
                double dis = fabs(rebuiltCoord[i*k + ki] - rebuiltCoord[j*k + ki]);
                chebyshev = dis>chebyshev ? dis : chebyshev;
            }
            chebyshevSum += chebyshev * 2;
        }
    }
    return chebyshevSum;
}

double SumDistance(const int k, const int n, const int dim, double* coord, int* pivots){
    //n:点数 k:支撑点个数
    // double* rebuiltCoord = (double*)malloc(sizeof(double) * n * k);
    // int i;
    // for(i=0; i<n*k; i++){
    //     rebuiltCoord[i] = 0;
    // }
    Distance(n, dim, coord);


    // Rebuild coordinates. New coordinate of one point is its distance to each pivot.
    // 计算基于支撑点的坐标的位置
    int i, ki;
    //对所有点进行遍历
    for(i=0; i<n; i++){
        //ki = 0, 1
        //新坐标的维数应该是k
        for(ki=0; ki<k; ki++){
            double distance = 0;
            //pivots[0] = 0  pivots[1] = i
            int pivoti = pivots[ki];

            rebuiltCoord[i*k + ki] = dist[pivoti * n + i];
            printf("\n%d\n", pivoti * n + i);
        }
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

void Circulation(const int k){
    int i, j;
    //初始化数组kk，使其从1开始递增
    for(i = k; i > 0; i --){
        for(int j = k; j > i; j --) kk[j] ++;
    }
}

void SortSum(const int k, const int n, const int dim, double* coord, int* pivots){

    double sum = SumDistance(k, n, dim, coord, pivots);
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < n; j ++) printf("%lf ", dist[i * n + j]);
        printf("\n");
    }

    for(int i = 0; i < n; i ++){
        for(int j = 0; j < k; j ++) printf("%lf ", rebuiltCoord[i * k + j]);
        printf("\n");
    }

    printf("%lf", sum);
}


void Combination(int ki, const int k, const int n, const int dim, const int M, double* coord, int* pivots,
                 double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots){
    if(ki==k-1){
        int i;

        for(i=pivots[ki-1]+1; i<n; i++){
            pivots[ki] = i;

            // Calculate sum of distance while combining different pivots.
            double distanceSum = SumDistance(k, n, dim, coord, pivots);

            // put data at the end of array
            maxDistanceSum[M] = distanceSum;
            minDistanceSum[M] = distanceSum;
            int kj;

            for(kj=0; kj<k; kj++){
                maxDisSumPivots[M*k + kj] = pivots[kj];
            }
            for(kj=0; kj<k; kj++){
                minDisSumPivots[M*k + kj] = pivots[kj];
            }

            // sort
            int a;
            
            for(a=M; a>0; a--){
                if(maxDistanceSum[a] > maxDistanceSum[a-1]){
                    double temp = maxDistanceSum[a];
                    maxDistanceSum[a] = maxDistanceSum[a-1];
                    maxDistanceSum[a-1] = temp;
                    int kj;
                    for(kj=0; kj<k; kj++){
                        int temp = maxDisSumPivots[a*k + kj];
                        maxDisSumPivots[a*k + kj] = maxDisSumPivots[(a-1)*k + kj];
                        maxDisSumPivots[(a-1)*k + kj] = temp;
                    }
                }
            }

              
            for(a=M; a>0; a--){
                if(minDistanceSum[a] < minDistanceSum[a-1]){
                    double temp = minDistanceSum[a];
                    minDistanceSum[a] = minDistanceSum[a-1];
                    minDistanceSum[a-1] = temp;
                    int kj;
                    for(kj=0; kj<k; kj++){
                        int temp = minDisSumPivots[a*k + kj];
                        minDisSumPivots[a*k + kj] = minDisSumPivots[(a-1)*k + kj];
                        minDisSumPivots[(a-1)*k + kj] = temp;
                    }
                }
            }
        }
        return;
    }

    // Recursively call Combination() to combine pivots
    int i; 
     

    for(i=pivots[ki-1]+1; i<n; i++) {
        pivots[ki] = i;
        
        Combination(ki+1, k, n, dim, M, coord, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);

        /** Iteration Log : pivots computed, best pivots, max distance sum, min distance sum pivots, min distance sum
        *** You can delete the logging code. **/
        if(ki==k-2){
            int kj;
            for(kj=0; kj<k; kj++){
                printf("%d ", pivots[kj]);
            }
            putchar('\t');
            for(kj=0; kj<k; kj++){
                printf("%d ", maxDisSumPivots[kj]);
            }
            printf("%lf\t", maxDistanceSum[0]);
            for(kj=0; kj<k; kj++){
                printf("%d ", minDisSumPivots[kj]);
            }
            printf("%lf\n", minDistanceSum[0]);
        }
    }
}

void main()
{
    int n = 5, dim = 2, k = 2;
    int pivots[2] = {0, 1};
    printf("dim = %d, n = %d, k = %d\n", dim, n, k);

    // Read Data
    double* coord = (double*)malloc(sizeof(double) * dim * n);
    int i;
    for(i=0; i<n; i++){
        int j;
        for(j=0; j<dim; j++){
            scanf("%lf", &coord[i*dim + j]);
        }
    }
    // Distance(n, dim, coord);
    kk = pow(n, k);
    SortSum(k, n, dim, coord, pivots);
    


}

